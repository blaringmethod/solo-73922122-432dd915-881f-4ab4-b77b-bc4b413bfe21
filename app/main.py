"""FastAPI application: versioned VCF upload + indexed interval queries.

Layout of the shared data directory:

    <data_dir>/meta.sqlite3          metadata, current-version pointer, ref pins
    <data_dir>/versions/<uuid>/      published file sets (variants.vcf.gz + .tbi)
    <data_dir>/staging/<instance>/<uuid>/   in-progress builds

Publishing is atomic: the file set is completed in staging, renamed into
versions/, and only then the current-version pointer is flipped inside one
SQLite transaction. Readers pin a version before opening its files, so an
in-flight query keeps its (version, data, index) trio consistent while a
replacement is published; the reaper frees old files once no pin remains.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import time
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass

import pysam
from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import RedirectResponse

from .db import Conflict, Store
from .lift import (
    FastaError,
    MappingError,
    check_mappings,
    lift_records,
    parse_fasta,
    parse_mappings,
)
from .vcfbuild import GZ_NAME, VCFError, validate_and_build

DEFAULT_MAX_UPLOAD = 10 * 1024 * 1024  # 10 MiB


def _fetch_intersecting(gz_path: str, chrom: str, start: int, end: int) -> list[str]:
    """Return raw records whose REF span intersects the 1-based closed interval.

    Fetching goes through the Tabix index; the per-record span check on top of
    the index iterator makes the closed-interval semantics exact (a deletion
    starting before `start` is still returned when it reaches into the range).
    """
    tb = pysam.TabixFile(gz_path)
    try:
        records = []
        for ln in tb.fetch(chrom, start - 1, end):
            f = ln.split("\t", 4)
            pos = int(f[1])
            if pos <= end and pos + len(f[3]) - 1 >= start:
                records.append(ln)
        return records
    finally:
        tb.close()


async def _read_limited(file: UploadFile, limit: int) -> bytes:
    """Read an upload in chunks, rejecting anything over `limit` bytes."""
    chunks = []
    size = 0
    while True:
        chunk = await file.read(1 << 20)
        if not chunk:
            break
        size += len(chunk)
        if size > limit:
            raise HTTPException(413, f"upload exceeds the {limit}-byte limit")
        chunks.append(chunk)
    return b"".join(chunks)


@dataclass
class Config:
    data_dir: str
    instance_id: str
    max_upload_bytes: int = DEFAULT_MAX_UPLOAD
    reaper_interval: float = 2.0      # seconds between reaper runs (0 disables)
    heartbeat_interval: float = 2.0   # seconds between heartbeats (0 disables)
    dead_after: float = 15.0          # instance considered dead after this silence
    orphan_grace: float = 60.0        # min age before orphan dirs are swept


def _config_from_env() -> Config:
    e = os.environ
    return Config(
        data_dir=e.get("VCF_DATA_DIR", "/data"),
        instance_id=e.get("VCF_INSTANCE_ID") or uuid.uuid4().hex,
        max_upload_bytes=int(e.get("VCF_MAX_UPLOAD_BYTES", DEFAULT_MAX_UPLOAD)),
        reaper_interval=float(e.get("VCF_REAPER_INTERVAL", "2")),
        heartbeat_interval=float(e.get("VCF_HEARTBEAT_INTERVAL", "2")),
        dead_after=float(e.get("VCF_DEAD_AFTER", "15")),
        orphan_grace=float(e.get("VCF_ORPHAN_GRACE", "60")),
    )


def create_app(config: Config | None = None) -> FastAPI:
    cfg = config or _config_from_env()
    store = Store(os.path.join(cfg.data_dir, "meta.sqlite3"))
    versions_dir = os.path.join(cfg.data_dir, "versions")
    staging_base = os.path.join(cfg.data_dir, "staging")
    staging_root = os.path.join(staging_base, cfg.instance_id)
    os.makedirs(versions_dir, exist_ok=True)
    # Leftover builds of a previous incarnation of this instance can never be
    # published any more; drop them before serving.
    shutil.rmtree(staging_root, ignore_errors=True)
    os.makedirs(staging_root, exist_ok=True)

    def reap_once() -> None:
        now = time.time()
        # 1) DB: drop dead instances' pins, delete reclaimable version rows.
        for _version, rel in store.reap_metadata(now, cfg.dead_after):
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        # 2) FS: sweep version dirs with no DB row (crashed publish/reaper).
        cutoff = now - cfg.orphan_grace
        known = store.all_version_paths()
        for name in os.listdir(versions_dir):
            path = os.path.join(versions_dir, name)
            rel = os.path.relpath(path, cfg.data_dir)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if rel not in known and old:
                shutil.rmtree(path, ignore_errors=True)
        # 3) FS: sweep staging dirs of dead instances.
        live = store.live_instances(now, cfg.dead_after)
        for inst in os.listdir(staging_base):
            path = os.path.join(staging_base, inst)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if inst not in live and old:
                shutil.rmtree(path, ignore_errors=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        store.heartbeat(cfg.instance_id, time.time())
        stopped = asyncio.Event()

        async def heartbeat_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.heartbeat_interval)
                try:
                    store.heartbeat(cfg.instance_id, time.time())
                except Exception:
                    pass

        async def reaper_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.reaper_interval)
                try:
                    await asyncio.to_thread(reap_once)
                except Exception:
                    pass

        tasks = []
        if cfg.heartbeat_interval > 0:
            tasks.append(asyncio.create_task(heartbeat_loop()))
        if cfg.reaper_interval > 0:
            tasks.append(asyncio.create_task(reaper_loop()))
        try:
            yield
        finally:
            stopped.set()
            for t in tasks:
                t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            try:
                store.unregister(cfg.instance_id)
            except Exception:
                pass

    app = FastAPI(
        title="Versioned VCF variant store",
        version="1.0.0",
        lifespan=lifespan,
        description=(
            "Upload single-chromosome biallelic VCFs (max 10 MiB, no sample "
            "columns, REF/ALT in ACGT); records are sorted, BGZF-compressed and "
            "Tabix-indexed. Query 1-based closed intervals; every record whose "
            "reference span [POS, POS+len(REF)-1] intersects the interval is "
            "returned via the index. Uploads are optimistic-concurrency "
            "controlled through `base_version`. `POST /lift` converts equal-"
            "length substitution records from the published assembly to a "
            "target reference through up to 100 gapless source->target "
            "interval mappings, honouring strand, and only reports a target "
            "variant when exactly one candidate matches the target FASTA."
        ),
    )
    app.state.cfg = cfg
    app.state.store = store
    app.state.reap_once = reap_once

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse("/docs")

    @app.post("/vcf", status_code=201)
    async def upload_vcf(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Version the client last saw (0 when empty)"),
    ):
        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = os.path.join(staging, "upload.vcf")
            size = 0
            with open(upload_path, "wb") as out:
                while True:
                    chunk = await file.read(1 << 20)
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > cfg.max_upload_bytes:
                        raise HTTPException(
                            413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                        )
                    out.write(chunk)
            try:
                chrom, n_records = await asyncio.to_thread(
                    validate_and_build, upload_path, staging
                )
            except VCFError as exc:
                raise HTTPException(400, f"invalid VCF: {exc}")
            os.remove(upload_path)

            # File set is complete; move it into place, then flip the pointer.
            final_rel = os.path.join("versions", uuid.uuid4().hex)
            os.rename(staging, os.path.join(cfg.data_dir, final_rel))
            try:
                new_version = await asyncio.to_thread(
                    store.publish, base_version, chrom, final_rel, n_records, time.time()
                )
            except Conflict as c:
                shutil.rmtree(os.path.join(cfg.data_dir, final_rel), ignore_errors=True)
                raise HTTPException(
                    409,
                    f"base_version {base_version} does not match current version "
                    f"{c.current}; re-read and retry",
                )
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
        return {"version": new_version, "chrom": chrom, "n_records": n_records}

    @app.get("/variants")
    def get_variants(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                records = _fetch_intersecting(gz_path, vchrom, start, end)
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "version": version,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin(ref_id)

    @app.post("/lift")
    async def lift(
        file: UploadFile = File(..., description="Single-chromosome target FASTA"),
        mappings: str = Form(
            ...,
            description=(
                "JSON array of at most 100 objects: "
                "{source_start, source_end, target_chrom, target_start, "
                "target_end, strand} with 1-based closed intervals, equal lengths, "
                "strand '+' or '-'."
            ),
        ),
        chrom: str = Query(..., description="Source chromosome as in the VCF"),
        start: int = Query(..., ge=1, description="1-based source interval start"),
        end: int = Query(..., ge=1, description="1-based source interval end"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")

        # --- validate the request body; any failure is a whole-request error ---
        try:
            parsed_mappings = parse_mappings(mappings)
        except MappingError as exc:
            raise HTTPException(400, f"invalid mappings: {exc}")

        raw = await _read_limited(file, cfg.max_upload_bytes)
        try:
            fasta_name, target_sequence = await asyncio.to_thread(
                parse_fasta, raw.decode("utf-8", errors="strict")
            )
        except UnicodeDecodeError as exc:
            raise HTTPException(400, f"invalid FASTA: not valid UTF-8 text: {exc}")
        except FastaError as exc:
            raise HTTPException(400, f"invalid FASTA: {exc}")
        try:
            check_mappings(
                parsed_mappings, fasta_name=fasta_name, target_length=len(target_sequence)
            )
        except MappingError as exc:
            raise HTTPException(400, f"invalid mappings: {exc}")

        # --- read the published data under a pin, so files cannot vanish ---
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                lines = await asyncio.to_thread(
                    _fetch_intersecting, gz_path, vchrom, start, end
                )
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            results = await asyncio.to_thread(
                lift_records, lines, parsed_mappings, target_sequence, version
            )
        finally:
            store.unpin(ref_id)

        lifted = [r for r in results if r["status"] == "lifted"]
        return {
            "version": version,
            "source": {
                "chrom": vchrom,
                "start": start,
                "end": end,
                "fasta_name": fasta_name,
            },
            "n_records": len(results),
            "n_lifted": len(lifted),
            "results": results,
        }

    @app.get("/version")
    def get_version():
        info = store.version_info()
        if info is None:
            return {"published": False, "current_version": 0}
        return {
            "published": True,
            "current_version": info["version"],
            "chrom": info["chrom"],
            "n_records": info["n_records"],
            "created_at": info["created_at"],
        }

    return app


# Run with: uvicorn --factory app.main:create_app
