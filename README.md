# Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement and safe concurrent access from two app instances.

## Run

```bash
docker compose up --build
```

Two instances come up sharing one data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published data; uploads may go to either one.

## API

### `POST /vcf?base_version=<n>` (multipart file field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T` (uppercase; no `N`, no
  symbolic alleles, no `*`, no multi-allelic `A,G`).

Records are sorted by position, BGZF-compressed and Tabix-indexed. The upload
carries the version the client last saw (`base_version`, `0` when the store is
empty). If two replacements race — even across the two instances — exactly one
wins; the loser gets `409` and must re-read `/version` and retry.
Returns `201` with `{"version", "chrom", "n_records"}`.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the
`version` that was actually read; that version's data and index are pinned for
the duration of the query, so version, data and index always match.

### `GET /version`

Current published version, chromosome and record count.

### `POST /lift?chrom=<c>&start=<s>&end=<e>` (multipart)

Converts substitution records from the published assembly to a new reference
("lifts" them), for a downstream annotation program.

Form fields:

- `file` — the target reference FASTA. At most 10 MiB, exactly one sequence,
  sequence bases limited to uppercase `A/C/G/T` (no `N`, no gaps, no lowercase).
  The sequence name is the target chromosome, `target_chrom` in the mappings
  must equal it, and every mapping must fit inside the sequence.
- `mappings` — JSON array of at most 100 objects:

  ```json
  [{"source_start": 1001, "source_end": 2000,
    "target_chrom": "chrT", "target_start": 5001, "target_end": 6000,
    "strand": "+"}]
  ```

  Coordinates are 1-based closed intervals. Each pair must have equal source
  and target lengths (gapless, no coordinate stretch). `strand` is `"+"` or
  `"-"`; on `"-"` the target interval is the reverse complement of the source
  interval. `target_chrom` must match the FASTA sequence name.

Query parameters are the 1-based closed *source* interval; records are read
through the existing Tabix index (validated VCF data only).

Only records whose REF and ALT have equal length (substitutions) are
converted. A record is converted when one single mapping fully covers its REF
span `[POS, POS+len(REF)-1]`; the converted target REF must then equal the
target FASTA at that position. A target variant is reported only when exactly
one *distinct* target placement passes that check (mappings that land on the
same target variant count as one). Every input record appears in `results`
with `status: "lifted"` and a `target`, or `status: "not_lifted"` plus a
`reason`:

| reason | meaning |
| --- | --- |
| `not_substitution` | REF and ALT have different lengths (indel/SV) |
| `no_coverage` | no mapping covers the REF span |
| `spans_segments` | the REF span crosses a mapping boundary |
| `reference_mismatch` | every candidate target placement disagrees with the FASTA |
| `multiple_targets` | more than one distinct target variant passed the check |

Each result entry carries the original raw `record`, the pinned
`vcf_version`, and the response carries the source interval. An invalid FASTA
or invalid mapping set fails the whole request with `400` and changes nothing;
oversize FASTA is `413`; no published data is `404`.

The lift reads through a version pin, so a concurrent replacement (from either
instance) neither changes the version/bases being lifted nor blocks the other
instance's upload and query.

```bash
curl -sS -X POST 'http://localhost:8001/lift?chrom=chr1&start=1&end=100000' \
  -F 'file=@target.fa' \
  -F 'mappings=[{"source_start":1,"source_end":1000,"target_chrom":"chrT",
       "target_start":5001,"target_end":6000,"strand":"+"}]'
```

## Consistency design

- **Atomic publish**: a new version is fully built (sorted, compressed,
  indexed) in a per-instance staging directory, moved into `versions/`, and
  only then the current-version pointer in SQLite is flipped inside a single
  `BEGIN IMMEDIATE` transaction. Readers never see a half-built file set.
- **Optimistic concurrency**: the pointer flip checks `base_version` inside
  the same transaction, so concurrent replacements from both instances are
  serialized and only one succeeds.
- **Old queries finish**: a query pins its version in a `refs` table before
  opening files. A background reaper deletes a superseded version's data and
  index only when no pin references it. Instances heartbeat; pins of crashed
  instances expire and never block reclamation.
- **Crash recovery**: if the process dies mid-build or mid-publish, the
  staging/orphan directories left behind are swept (with a grace period) and
  the previously published version keeps serving. The pointer flip is the last
  step, so a restart never exposes a partial file set.

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 90 tests
uvicorn --factory app.main:create_app   # single dev instance (uses /data or $VCF_DATA_DIR)
```

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`, `VCF_HEARTBEAT_INTERVAL`,
`VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
