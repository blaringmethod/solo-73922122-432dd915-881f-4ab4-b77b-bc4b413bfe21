"""Tests for POST /lift: coordinate conversion, strand handling, reasons."""
import json
import os

import pytest
from conftest import SAMPLE, upload

import app.main as main_module
from app.lift import parse_fasta

FASTA_NAME = "chrT"


def make_fasta(length=2000, overlays=None, name=FASTA_NAME, line_width=None):
    """Build a FASTA string; `overlays` maps 1-based pos -> base string."""
    seq = ["A"] * length
    for pos, bases in (overlays or {}).items():
        for i, base in enumerate(bases):
            seq[pos - 1 + i] = base
    text = "".join(seq)
    header = f">{name}\n"
    if line_width:
        body = "\n".join(
            text[i : i + line_width] for i in range(0, len(text), line_width)
        )
    else:
        body = text
    return header + body + "\n"


def mapping(source_start, source_end, target_start, target_end, strand="+",
            target_chrom=FASTA_NAME):
    return {
        "source_start": source_start,
        "source_end": source_end,
        "target_chrom": target_chrom,
        "target_start": target_start,
        "target_end": target_end,
        "strand": strand,
    }


def lift(client, fasta, mappings, chrom="chr1", start=1, end=1000000):
    body = mappings if isinstance(mappings, str) else json.dumps(mappings)
    return client.post(
        "/lift",
        params={"chrom": chrom, "start": start, "end": end},
        files={"file": ("target.fa", fasta.encode() if isinstance(fasta, str) else fasta)},
        data={"mappings": body},
    )


def entry_for(body, pos):
    matches = [r for r in body["results"] if r["pos"] == pos]
    assert len(matches) == 1, f"expected one entry for source pos {pos}"
    return matches[0]


SRC = (
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tG\t.\t.\t.\n"
    "chr1\t150\t.\tACG\tTTG\t.\t.\t.\n"
    "chr1\t300\t.\tT\tC\t.\t.\t.\n"
)


@pytest.fixture()
def src_client(client):
    r = upload(client, SRC, base_version=0)
    assert r.status_code == 201, r.text
    return client


# ------------------------------------------------------------------ happy paths


def test_plus_strand_lifts_each_substitution(src_client):
    fasta = make_fasta(overlays={1050: "A", 1100: "ACG", 1250: "T"})
    r = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["version"] == 1
    assert body["n_records"] == 3 and body["n_lifted"] == 3

    e = entry_for(body, 100)
    assert e["status"] == "lifted"
    assert e["target"] == {
        "chrom": "chrT", "pos": 1050, "ref": "A", "alt": "G",
        "strand": "+", "mapping_index": 0,
    }

    e = entry_for(body, 150)
    assert e["target"]["pos"] == 1100
    assert (e["target"]["ref"], e["target"]["alt"]) == ("ACG", "TTG")

    e = entry_for(body, 300)
    assert e["target"] == {
        "chrom": "chrT", "pos": 1250, "ref": "T", "alt": "C",
        "strand": "+", "mapping_index": 0,
    }


def test_minus_strand_reverse_complements_coordinates_and_bases(src_client):
    # length of mapping 301; target offset = 301 - source offset - len(REF)
    fasta = make_fasta(overlays={750: "T", 698: "CGT", 550: "A"})
    r = lift(src_client, fasta, [mapping(50, 350, 500, 800, strand="-")])
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_lifted"] == 3

    e = entry_for(body, 100)  # len 1
    assert e["target"]["pos"] == 750
    assert (e["target"]["ref"], e["target"]["alt"]) == ("T", "C")

    e = entry_for(body, 150)  # ACG -> revcomp CGT ; TTG -> CAA
    assert e["target"]["pos"] == 698
    assert (e["target"]["ref"], e["target"]["alt"]) == ("CGT", "CAA")

    e = entry_for(body, 300)
    assert e["target"]["pos"] == 550
    assert (e["target"]["ref"], e["target"]["alt"]) == ("A", "G")


def test_mapping_boundaries_are_inclusive(src_client):
    """A record at exactly source_start / source_end is covered."""
    body_vcf = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t50\t.\tA\tG\t.\t.\t.\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\n"
    )
    upload(src_client, body_vcf, base_version=1)
    fasta = make_fasta(overlays={1000: "A", 1050: "A"})
    r = lift(src_client, fasta, [mapping(50, 100, 1000, 1050)])
    body = r.json()
    assert body["n_lifted"] == 2
    assert entry_for(body, 50)["target"]["pos"] == 1000
    assert entry_for(body, 100)["target"]["pos"] == 1050


def test_each_result_carries_original_record_and_pinned_version(src_client):
    fasta = make_fasta(overlays={1050: "A"})
    body = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)]).json()
    e = entry_for(body, 100)
    assert e["record"] == "chr1\t100\t.\tA\tG\t.\t.\t."
    assert e["vcf_version"] == 1
    assert body["version"] == 1
    assert body["source"]["chrom"] == "chr1"
    assert body["source"]["fasta_name"] == FASTA_NAME


def test_query_interval_limits_records_via_index(src_client):
    fasta = make_fasta(overlays={1050: "A", 1100: "ACG", 1250: "T"})
    r = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)], start=120, end=160)
    body = r.json()
    assert [x["pos"] for x in body["results"]] == [150]


def test_multiline_fasta_and_description_in_header(src_client):
    fasta = make_fasta(overlays={1050: "A"}, line_width=7)
    fasta = fasta.replace(">chrT", ">chrT assembled target", 1)
    r = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 200, r.text
    assert entry_for(r.json(), 100)["status"] == "lifted"


def test_crlf_fasta_is_accepted(src_client):
    fasta = make_fasta(overlays={1050: "A"}).replace("\n", "\r\n")
    r = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 200, r.text


# ------------------------------------------------------------- per-record reasons


def test_insertion_and_deletion_are_not_substitutions(client):
    upload(client, SAMPLE, base_version=0)
    r = lift(client, make_fasta(), [mapping(1, 1000, 1, 1000)])
    body = r.json()
    for pos in (100, 110, 200):
        e = entry_for(body, pos)
        assert e["status"] == "not_lifted"
        assert e["reason"] == "not_substitution"
        assert e["target"] is None
    # insertion at 110 is REF=C ALT=CT -> unequal lengths
    assert entry_for(body, 110)["reason"] == "not_substitution"


def test_no_coverage_reports_reason(src_client):
    r = lift(src_client, make_fasta(), [mapping(200, 400, 1000, 1200)])
    e = entry_for(r.json(), 100)
    assert e["status"] == "not_lifted"
    assert e["reason"] == "no_coverage"
    assert e["target"] is None
    assert "100" in e["detail"]


def test_record_spanning_a_segment_boundary_reports_spans_segments(src_client):
    # two adjacent source segments 50-149 and 150-250; REF at 148 len 3 -> 148-150
    upload(
        src_client,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t148\t.\tACG\tTTG\t.\t.\t.\n",
        base_version=1,
    )
    fasta = make_fasta()
    r = lift(
        src_client,
        fasta,
        [mapping(50, 149, 1000, 1099), mapping(150, 250, 1200, 1300)],
    )
    body = r.json()
    e = entry_for(body, 148)
    assert e["reason"] == "spans_segments"
    assert e["target"] is None
    assert "148-150" in e["detail"]


def test_all_candidates_fail_reference_check(src_client):
    # mapping points at a target base that does not match REF
    fasta = make_fasta(overlays={1050: "C"})
    body = lift(src_client, fasta, [mapping(50, 350, 1000, 1300)]).json()
    e = entry_for(body, 100)
    assert e["status"] == "not_lifted"
    assert e["reason"] == "reference_mismatch"
    assert e["target"] is None


def test_multiple_distinct_targets_both_matching_is_refused(src_client):
    fasta = make_fasta(overlays={150: "A", 550: "A"})
    # two different source segments (same record fully inside each) -> different targets
    r = lift(
        src_client,
        fasta,
        [mapping(50, 350, 100, 400), mapping(50, 350, 500, 800)],
    )
    body = r.json()
    e = entry_for(body, 100)
    assert e["status"] == "not_lifted"
    assert e["reason"] == "multiple_targets"
    assert e["target"] is None
    assert "2" in e["detail"]


def test_one_matching_target_wins_when_other_candidates_mismatch(src_client):
    # mapping 0 lands at target 1050 where FASTA has C; mapping 1 at 1250 has A
    fasta = make_fasta(overlays={1050: "C", 1300: "A"})
    body = lift(
        src_client,
        fasta,
        [mapping(50, 350, 1000, 1300), mapping(50, 350, 1250, 1550)],
    ).json()
    e = entry_for(body, 100)
    assert e["status"] == "lifted"
    assert e["target"]["mapping_index"] == 1
    assert e["target"]["pos"] == 1300


def test_duplicate_mappings_to_same_target_are_one_variant(src_client):
    fasta = make_fasta(overlays={1050: "A"})
    body = lift(
        src_client,
        fasta,
        [mapping(50, 350, 1000, 1300), mapping(50, 350, 1000, 1300)],
    ).json()
    e = entry_for(body, 100)
    assert e["status"] == "lifted", "identical target placements are not ambiguous"
    assert e["target"]["pos"] == 1050


# --------------------------------------------------------------- request errors


BAD_FASTAS = {
    "missing_header": "ACGTACGT\n",
    "empty_sequence": ">chrT\n",
    "header_only_blank": ">chrT\n\n",
    "two_sequences": ">chrT\nACGT\n>chrU\nACGT\n",
    "lowercase": ">chrT\nacgt\n",
    "with_n": ">chrT\nACGTN\n",
    "with_gap": ">chrT\nACGT-\n",
    "with_digit": ">chrT\nACGT1\n",
    "empty_name": ">\nACGT\n",
    "data_before_header": "ACGT\n>chrT\nACGT\n",
    "sequence_with_space": ">chrT\nAC GT\n",
}


@pytest.mark.parametrize("name", sorted(BAD_FASTAS))
def test_invalid_fasta_rejects_whole_request(src_client, app, name):
    before = app.state.store.version_info()
    r = lift(src_client, BAD_FASTAS[name], [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 400, f"{name}: {r.text}"
    assert "invalid FASTA" in r.json()["detail"]
    assert app.state.store.version_info() == before
    q = src_client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
    assert q.status_code == 200 and q.json()["version"] == 1


BAD_MAPPINGS = {
    "not_json": "not json",
    "not_array": '{"source_start": 1}',
    "empty_object": [{}],
    "unknown_field": [dict(mapping(50, 350, 1000, 1300), extra=1)],
    "unequal_lengths": [mapping(50, 350, 1000, 1301)],
    "reversed_source": [mapping(350, 50, 1000, 1300)],
    "reversed_target": [mapping(50, 350, 1300, 1000)],
    "bad_strand": [dict(mapping(50, 350, 1000, 1300), strand="*")],
    "strand_missing": [{k: v for k, v in mapping(50, 350, 1000, 1300).items() if k != "strand"}],
    "float_coord": [dict(mapping(50, 350, 1000, 1300), source_start=50.5)],
    "string_coord": [dict(mapping(50, 350, 1000, 1300), source_start="50")],
    "bool_coord": [dict(mapping(50, 350, 1000, 1300), source_start=True)],
    "zero_coord": [mapping(0, 350, 1000, 1300)],
    "target_chrom_mismatch": [mapping(50, 350, 1000, 1300, target_chrom="chrU")],
    "target_beyond_fasta": [mapping(50, 350, 1900, 2200)],
}


@pytest.mark.parametrize("name", sorted(BAD_MAPPINGS))
def test_invalid_mappings_reject_whole_request(src_client, app, name):
    before = app.state.store.version_info()
    r = lift(src_client, make_fasta(), BAD_MAPPINGS[name])
    assert r.status_code == 400, f"{name}: {r.text}"
    assert "invalid mappings" in r.json()["detail"]
    assert app.state.store.version_info() == before
    assert src_client.get(
        "/variants", params={"chrom": "chr1", "start": 1, "end": 1000}
    ).status_code == 200


def test_more_than_100_mappings_rejected(src_client):
    maps = [mapping(50, 350, 1000, 1300) for _ in range(101)]
    r = lift(src_client, make_fasta(), maps)
    assert r.status_code == 400
    assert "at most 100" in r.json()["detail"]


def test_exactly_100_mappings_accepted(src_client):
    maps = [mapping(50, 350, 1000, 1300) for _ in range(100)]
    fasta = make_fasta(overlays={1050: "A"})
    assert lift(src_client, fasta, maps).status_code == 200


def test_empty_mapping_list_is_valid_but_lifts_nothing(src_client):
    body = lift(src_client, make_fasta(), []).json()
    assert body["n_records"] == 3 and body["n_lifted"] == 0
    assert {r["reason"] for r in body["results"]} == {"no_coverage"}


def test_oversize_fasta_rejected(client, app):
    big = ">chrT\n" + "A" * (app.state.cfg.max_upload_bytes + 10) + "\n"
    upload(client, SRC, base_version=0)
    r = lift(client, big, [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 413
    assert client.get("/version").json()["current_version"] == 1


def test_non_utf8_fasta_rejected(src_client):
    r = lift(src_client, b">chrT\n\xff\xfe\n", [mapping(50, 350, 1000, 1300)])
    assert r.status_code == 400
    assert "invalid FASTA" in r.json()["detail"]


# ------------------------------------------------------------- state / lifecycle


def test_lift_without_published_data_is_404(client):
    r = lift(client, make_fasta(), [mapping(1, 10, 1, 10)])
    assert r.status_code == 404


def test_lift_wrong_chrom_is_400(src_client):
    r = lift(src_client, make_fasta(), [mapping(50, 350, 1000, 1300)], chrom="chr2")
    assert r.status_code == 400
    assert "chr1" in r.json()["detail"]


def test_lift_reversed_interval_is_400(src_client):
    r = lift(src_client, make_fasta(), [mapping(50, 350, 1000, 1300)], start=500, end=100)
    assert r.status_code == 400


def test_lift_outside_query_interval_has_no_records(src_client):
    body = lift(
        src_client, make_fasta(), [mapping(50, 350, 1000, 1300)], start=500, end=600
    ).json()
    assert body["n_records"] == 0 and body["results"] == []


def test_lift_uses_pinned_version_during_concurrent_replacement(tmp_path):
    """A replacement mid-lift must not change the version/bases being lifted,
    and the other instance's upload and query must stay unblocked."""
    from fastapi.testclient import TestClient

    from app.main import Config, create_app

    app1 = create_app(Config(data_dir=str(tmp_path / "d"), instance_id="a",
                             reaper_interval=0, heartbeat_interval=0))
    app2 = create_app(Config(data_dir=str(tmp_path / "d"), instance_id="b",
                             reaper_interval=0, heartbeat_interval=0))
    with TestClient(app1) as c1, TestClient(app2) as c2:
        assert upload(c1, SRC, base_version=0).status_code == 201

        original = main_module._fetch_intersecting
        seen = {}
        calls = {"n": 0}

        def fetch_then_publish(*args, **kwargs):
            # Publish version 2 from the *other* instance while lift is running.
            # Guarded so the re-entrant query/variants calls do not recurse.
            calls["n"] += 1
            if calls["n"] == 1:
                up = upload(c2, "##fileformat=VCFv4.2\n"
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                                "chr1\t900\t.\tA\tG\t.\t.\t.\n", base_version=1)
                seen["upload"] = up.status_code
                seen["query"] = c2.get(
                    "/variants", params={"chrom": "chr1", "start": 1, "end": 2000}
                ).status_code
            return original(*args, **kwargs)

        main_module._fetch_intersecting = fetch_then_publish
        try:
            fasta = make_fasta(overlays={1050: "A"})
            r = lift(c1, fasta, [mapping(50, 350, 1000, 1300)])
        finally:
            main_module._fetch_intersecting = original

        assert seen == {"upload": 201, "query": 200}
        body = r.json()
        # lift still reports the version it pinned, and the original records
        assert body["version"] == 1
        assert [x["pos"] for x in body["results"]] == [100, 150, 300]
        assert entry_for(body, 100)["vcf_version"] == 1
        # the new version is visible afterwards
        assert c1.get("/version").json()["current_version"] == 2


def test_published_data_unchanged_after_invalid_lift(src_client, app):
    info_before = app.state.store.version_info()
    data_dir = os.path.join(app.state.cfg.data_dir, info_before["path"])
    files_before = sorted(os.listdir(data_dir))
    lift(src_client, ">chrT\nACGTN\n", [mapping(50, 350, 1000, 1300)])
    assert app.state.store.version_info() == info_before
    assert sorted(os.listdir(data_dir)) == files_before


# ------------------------------------------------------------------- docs / meta


def test_openapi_documents_lift(client):
    spec = client.get("/openapi.json").json()
    assert "/lift" in spec["paths"]
    assert client.get("/docs").status_code == 200


def test_parse_fasta_direct():
    name, seq = parse_fasta(">chr1 desc\nACGT\nACGT\n")
    assert name == "chr1"
    assert seq == "ACGTACGT"
