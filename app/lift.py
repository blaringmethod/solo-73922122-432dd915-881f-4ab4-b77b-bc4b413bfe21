"""Variant lifting from a source assembly to a target reference.

This module is pure: it validates the target FASTA and the user supplied
source->target mappings, then converts VCF substitution records and checks
each candidate target reference against the target FASTA.  File I/O, VCF
index reads and version pinning live in :mod:`app.main`.

A mapping is a *gapless, equal-length* closed interval pair::

    source [source_start, source_end]  ->  target [target_start, target_end]
    with source_end - source_start == target_end - target_start

``strand`` is ``"+"`` (same orientation) or ``"-"`` (target interval is the
reverse complement of the source interval).  The source chromosome is the
published VCF's chromosome; the target chromosome is named explicitly because
the target FASTA is a separate reference.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass

MAX_MAPPINGS = 100

# The target FASTA is deliberately strict: uppercase A/C/G/T only.
ACGT_RE = re.compile(r"[ACGT]+")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Machine-readable reasons for records that are not lifted.
REASON_NOT_SUBSTITUTION = "not_substitution"
REASON_NO_COVERAGE = "no_coverage"
REASON_SPANS_SEGMENTS = "spans_segments"
REASON_REFERENCE_MISMATCH = "reference_mismatch"
REASON_MULTIPLE_TARGETS = "multiple_targets"

STATUS_LIFTED = "lifted"
STATUS_NOT_LIFTED = "not_lifted"

_REQUIRED_MAPPING_KEYS = (
    "source_start",
    "source_end",
    "target_chrom",
    "target_start",
    "target_end",
    "strand",
)


class FastaError(ValueError):
    """The uploaded target FASTA violates the acceptance rules."""


class MappingError(ValueError):
    """The supplied mapping set is malformed."""


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase A/C/G/T sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def parse_fasta(text: str) -> tuple[str, str]:
    """Validate a single-chromosome FASTA and return ``(name, sequence)``.

    Rules: exactly one ``>`` header, a non-empty sequence name, at least one
    sequence base, and only uppercase ``A/C/G/T`` in the sequence.  Blank
    lines are ignored; anything else is rejected.
    """
    name: str | None = None
    parts: list[str] = []
    saw_header = False

    for lineno, line in enumerate(text.splitlines(), 1):
        if line == "":
            continue
        if line.startswith(">"):
            if saw_header:
                raise FastaError(
                    f"line {lineno}: multiple sequences found; the target FASTA "
                    "must contain exactly one chromosome"
                )
            header = line[1:].strip()
            if not header:
                raise FastaError(f"line {lineno}: sequence name is empty")
            name = header.split()[0]
            saw_header = True
            continue
        if not saw_header:
            raise FastaError(f"line {lineno}: sequence data before the '>' header")
        if not ACGT_RE.fullmatch(line):
            raise FastaError(
                f"line {lineno}: sequence must contain only uppercase A/C/G/T"
            )
        parts.append(line)

    if not saw_header:
        raise FastaError("missing '>' header line")
    sequence = "".join(parts)
    if not sequence:
        raise FastaError("sequence is empty")
    assert name is not None  # implied by saw_header
    return name, sequence


@dataclass(frozen=True)
class Mapping:
    index: int
    source_start: int
    source_end: int
    target_chrom: str
    target_start: int
    target_end: int
    strand: str

    @property
    def length(self) -> int:
        return self.source_end - self.source_start + 1


@dataclass(frozen=True)
class Placement:
    mapping_index: int
    chrom: str
    pos: int
    ref: str
    alt: str
    strand: str


def _parse_int(obj: dict, key: str, prefix: str) -> int:
    value = obj[key]
    # bool is an int subclass; coordinates must be actual JSON integers.
    if type(value) is not int:
        raise MappingError(f"{prefix}: '{key}' must be an integer")
    if value < 1:
        raise MappingError(f"{prefix}: '{key}' must be >= 1")
    return value


def _parse_mapping(index: int, obj: object) -> Mapping:
    prefix = f"mapping {index}"
    if not isinstance(obj, dict):
        raise MappingError(f"{prefix}: must be a JSON object")

    unknown = sorted(set(obj) - set(_REQUIRED_MAPPING_KEYS))
    if unknown:
        raise MappingError(f"{prefix}: unknown field(s): {', '.join(unknown)}")
    missing = [k for k in _REQUIRED_MAPPING_KEYS if k not in obj]
    if missing:
        raise MappingError(f"{prefix}: missing field(s): {', '.join(missing)}")

    source_start = _parse_int(obj, "source_start", prefix)
    source_end = _parse_int(obj, "source_end", prefix)
    target_start = _parse_int(obj, "target_start", prefix)
    target_end = _parse_int(obj, "target_end", prefix)

    if source_start > source_end:
        raise MappingError(f"{prefix}: source_start must be <= source_end")
    if target_start > target_end:
        raise MappingError(f"{prefix}: target_start must be <= target_end")
    source_len = source_end - source_start + 1
    target_len = target_end - target_start + 1
    if source_len != target_len:
        raise MappingError(
            f"{prefix}: source interval length {source_len} != target interval "
            f"length {target_len}; mappings must be equal-length and gapless"
        )

    target_chrom = obj["target_chrom"]
    if not isinstance(target_chrom, str) or not target_chrom:
        raise MappingError(f"{prefix}: 'target_chrom' must be a non-empty string")

    strand = obj["strand"]
    if strand not in ("+", "-"):
        raise MappingError(f"{prefix}: 'strand' must be '+' or '-'")

    return Mapping(
        index=index,
        source_start=source_start,
        source_end=source_end,
        target_chrom=target_chrom,
        target_start=target_start,
        target_end=target_end,
        strand=strand,
    )


def parse_mappings(raw: str) -> list[Mapping]:
    """Parse and structurally validate the JSON mapping array."""
    try:
        data = json.loads(raw)
    except (TypeError, json.JSONDecodeError) as exc:
        raise MappingError(f"must be a JSON array of mapping objects: {exc}") from exc
    if not isinstance(data, list):
        raise MappingError("must be a JSON array of mapping objects")
    if len(data) > MAX_MAPPINGS:
        raise MappingError(
            f"at most {MAX_MAPPINGS} mappings are allowed, got {len(data)}"
        )
    return [_parse_mapping(i, item) for i, item in enumerate(data)]


def check_mappings(
    mappings: list[Mapping], *, fasta_name: str, target_length: int
) -> None:
    """Validate mappings against the target FASTA (name and coordinate bounds)."""
    for m in mappings:
        if m.target_chrom != fasta_name:
            raise MappingError(
                f"mapping {m.index}: target_chrom '{m.target_chrom}' does not "
                f"match the target FASTA sequence '{fasta_name}'"
            )
        if m.target_end > target_length:
            raise MappingError(
                f"mapping {m.index}: target interval ends at {m.target_end}, "
                f"beyond the target sequence length {target_length}"
            )


def _placement(mapping: Mapping, pos: int, ref: str, alt: str) -> Placement:
    """Convert one fully covered source span to a target placement."""
    offset = pos - mapping.source_start  # 0-based offset into the source interval
    if mapping.strand == "+":
        target_pos = mapping.target_start + offset
        target_ref, target_alt = ref, alt
    else:
        # Reverse orientation: source offset o maps to target offset
        # (length - o - len(ref)) .. (length - o - 1).
        target_pos = mapping.target_start + mapping.length - offset - len(ref)
        target_ref, target_alt = revcomp(ref), revcomp(alt)
    return Placement(
        mapping_index=mapping.index,
        chrom=mapping.target_chrom,
        pos=target_pos,
        ref=target_ref,
        alt=target_alt,
        strand=mapping.strand,
    )


def lift_record(
    line: str,
    mappings: list[Mapping],
    target_sequence: str,
    source_version: int,
) -> dict:
    """Lift one raw VCF record, returning a JSON-ready result entry.

    The target variant is only reported when exactly one *distinct* candidate
    placement passes the target-reference check.  `source_version` is the
    pinned VCF version the record was read from; it is copied verbatim onto
    every entry so a concurrent publish cannot desynchronise the citation.
    """
    fields = line.split("\t")
    if len(fields) < 5:  # pragma: no cover - published records are validated
        raise ValueError("published VCF record has fewer than 5 columns")
    chrom, pos_s, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
    pos = int(pos_s)

    result = {
        "record": line,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "vcf_version": source_version,
        "status": STATUS_NOT_LIFTED,
        "reason": None,
        "detail": None,
        "target": None,
    }

    # Only substitutions (REF and ALT of equal length) carry over by simple
    # coordinate transformation; indels/SVs are reported, never guessed.
    if len(ref) != len(alt):
        result["reason"] = REASON_NOT_SUBSTITUTION
        result["detail"] = (
            f"REF length {len(ref)} != ALT length {len(alt)}; only equal-length "
            "substitution records are lifted"
        )
        return result

    span_end = pos + len(ref) - 1
    containing = [m for m in mappings if m.source_start <= pos and m.source_end >= span_end]
    if not containing:
        overlapping = [
            m for m in mappings if m.source_start <= span_end and m.source_end >= pos
        ]
        if overlapping:
            result["reason"] = REASON_SPANS_SEGMENTS
            result["detail"] = (
                f"REF span {pos}-{span_end} is not fully contained in any single "
                f"mapping (partially covered by mapping(s) "
                f"{', '.join(str(m.index) for m in overlapping)})"
            )
        else:
            result["reason"] = REASON_NO_COVERAGE
            result["detail"] = f"REF span {pos}-{span_end} is not covered by any mapping"
        return result

    # Distinct target variants: overlapping mappings that yield the exact same
    # target variant are one placement, not an ambiguity.
    candidates: dict[tuple, Placement] = {}
    for m in containing:
        p = _placement(m, pos, ref, alt)
        candidates.setdefault((p.chrom, p.pos, p.ref, p.alt), p)

    matched: list[Placement] = []
    for p in candidates.values():
        observed = target_sequence[p.pos - 1 : p.pos - 1 + len(p.ref)]
        if observed == p.ref:
            matched.append(p)

    if len(matched) == 1:
        p = matched[0]
        result["status"] = STATUS_LIFTED
        result["detail"] = (
            f"mapped by mapping {p.mapping_index} ({p.strand} strand)"
        )
        result["target"] = {
            "chrom": p.chrom,
            "pos": p.pos,
            "ref": p.ref,
            "alt": p.alt,
            "strand": p.strand,
            "mapping_index": p.mapping_index,
        }
    elif not matched:
        result["reason"] = REASON_REFERENCE_MISMATCH
        result["detail"] = (
            f"{len(candidates)} candidate target placement(s) checked; the target "
            "FASTA reference did not match any of them"
        )
    else:
        result["reason"] = REASON_MULTIPLE_TARGETS
        result["detail"] = (
            f"{len(matched)} distinct target variants passed the reference check; "
            "refusing to choose one"
        )
    return result


def lift_records(
    lines: list[str],
    mappings: list[Mapping],
    target_sequence: str,
    source_version: int,
) -> list[dict]:
    """Lift a batch of raw VCF records, preserving input order."""
    return [lift_record(line, mappings, target_sequence, source_version) for line in lines]
