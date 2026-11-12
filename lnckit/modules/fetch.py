"""SRA fetch stage - accession list -> local FASTQ, the true start of the pipeline.

Biological purpose
------------------
Re-analyses of published plant RNA-seq start from run accessions (SRR/ERR/DRR).
Making this a pipeline stage (rather than a manual `fastq-dump` loop) means the
accessions, the split mode and the toolkit version all land in `run_metadata.json`
alongside everything else, so the run is reproducible from the accession list up.

Inputs   : `inputs.sra_accessions` (list) and/or `inputs.sra_accession_file` (one
           accession per line, `#` comments allowed).
Outputs  : gzipped FASTQ in `<output>/raw_reads/`, named `<ACC>_1.fastq.gz` /
           `<ACC>_2.fastq.gz` (paired) or `<ACC>.fastq.gz` (single).
Tools    : sra-tools (`prefetch` + `fastq-dump` or `fasterq-dump`).
Runtime  : network-bound; typically minutes per run, tens of minutes for deep runs.
Resume   : idempotent - an accession whose FASTQ already exists is skipped.
"""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import List, Optional

from lnckit.core.context import RunContext
from lnckit.core.exceptions import OutputValidationError
from lnckit.core.logging_utils import get_logger
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.fetch")

ACCESSION_PREFIXES = ("SRR", "ERR", "DRR", "SRX", "ERX", "DRX")


def read_accessions(cfg) -> List[str]:
    """Collect accessions from the inline list and/or an accession file."""
    accs: List[str] = list(cfg.get("inputs.sra_accessions") or [])
    acc_file = cfg.get("inputs.sra_accession_file")
    if acc_file:
        for line in Path(acc_file).read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                accs.append(line)
    seen, unique = set(), []
    for a in accs:
        a = a.strip()
        if a and a not in seen:
            seen.add(a)
            unique.append(a)
    for a in unique:
        if not a.upper().startswith(ACCESSION_PREFIXES):
            log.warning("'%s' does not look like an SRA run accession "
                        "(expected SRR/ERR/DRR...)", a)
    return unique


def existing_fastqs(reads_dir: Path, acc: str) -> List[Path]:
    """FASTQ files already present for an accession (paired, single or unpaired)."""
    return sorted(p for p in reads_dir.glob(f"{acc}*.fastq*")
                  if not p.name.endswith(".tmp"))


def _fetch_marker(reads_dir: Path, acc: str) -> Path:
    """Sentinel written once an accession's FASTQ(s) are known-complete."""
    return reads_dir / f".{acc}.fetched"


def _gzip_ok(path: Path) -> bool:
    """True if a .gz decompresses cleanly to EOF - catches truncated/partial downloads
    (the ones fastp later rejects with 'igzip: unexpected eof')."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def _existing_is_complete(reads_dir: Path, acc: str, already: List[Path]) -> bool:
    """Decide whether already-present FASTQs may be reused (resume) or must be re-fetched.

    Fast path: a `.fetched` marker means we verified them before. Otherwise verify each
    gzip once and, if sound, drop the marker so later resumes are instant; a truncated
    file (from an interrupted earlier fetch) is reported so the caller re-downloads.
    """
    marker = _fetch_marker(reads_dir, acc)
    if marker.exists():
        return True
    bad = [p for p in already if p.name.endswith(".gz") and not _gzip_ok(p)]
    if bad:
        log.warning("%s: incomplete/corrupt FASTQ from an earlier run (%s) - re-downloading",
                    acc, ", ".join(p.name for p in bad))
        return False
    marker.touch()                        # remember it's good -> O(1) on future resumes
    return True


def _gzip_in_place(path: Path) -> Path:
    """Compress an uncompressed FASTQ (fasterq-dump cannot write .gz itself)."""
    target = path.with_suffix(path.suffix + ".gz")
    with open(path, "rb") as fin, gzip.open(target, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    path.unlink()
    return target


def _prefetched_sra(cache_dir: Path, acc: str) -> Optional[Path]:
    """Locate the .sra archive prefetch wrote (layout differs across toolkit versions)."""
    for candidate in (cache_dir / acc / f"{acc}.sra", cache_dir / f"{acc}.sra"):
        if candidate.exists():
            return candidate
    hits = list(cache_dir.glob(f"**/{acc}.sra"))
    return hits[0] if hits else None


def _clear_stale_locks(cache_dir: Path, acc: str) -> None:
    """Remove leftover prefetch `.lock` files from an interrupted earlier download.

    prefetch writes `<acc>.sra.lock` while downloading and deletes it on success; a
    crash/kill/disconnect leaves it behind, and the next prefetch then refuses ("lock
    exists ... download canceled"). LncKit fetches ONE accession at a time in a single
    process, so any lock we find here is stale and safe to remove before retrying.
    """
    for lock in cache_dir.glob(f"**/{acc}*.lock"):
        try:
            lock.unlink()
            log.info("removed stale prefetch lock: %s", lock)
        except OSError:
            pass


def _purge_sra_cache(cache_dir: Path, acc: str) -> None:
    """Delete an accession's cached .sra - its per-accession dir plus any flat-layout or
    partial/.tmp/.lock siblings - so an interrupted download is re-fetched from scratch
    instead of being dumped half-complete."""
    d = cache_dir / acc
    if d.is_dir():
        shutil.rmtree(d, ignore_errors=True)
    for f in cache_dir.glob(f"{acc}.sra*"):
        try:
            f.unlink()
        except OSError:
            pass


def run_fetch(ctx: RunContext) -> Path:
    """Download every configured accession into `<output>/raw_reads/`.

    Returns the reads directory, which becomes the pipeline's input directory.
    """
    cfg = ctx.config
    accessions = read_accessions(cfg)
    # path only - raw_reads/ is created just before the first dump (below), so during the
    # download step only sra_cache/ appears; raw_reads/ shows up when fastq-dump starts.
    reads_dir = ctx.output_dir / "raw_reads"

    if not accessions:
        log.info("no SRA accessions configured - using inputs.reads_dir as-is")
        return Path(cfg.get("inputs.reads_dir") or reads_dir)

    sra_cfg = cfg.get("sra", {}) or {}
    dumper = sra_cfg.get("tool", "fastq-dump")
    split = sra_cfg.get("split", "split-files")
    use_prefetch = sra_cfg.get("prefetch", True)
    want_gzip = sra_cfg.get("gzip", True)
    cache_dir = ctx.dirs["sra_cache"]

    log.info("fetching %d accession(s) with %s", len(accessions), dumper)

    for acc in accessions:
        already = existing_fastqs(reads_dir, acc)
        if already and cfg.get("execution.resume", True) \
                and _existing_is_complete(reads_dir, acc, already):
            log.info("%s: %d FASTQ file(s) already present - skipping",
                     acc, len(already))
            continue
        if already:                       # stale/corrupt from an interrupted run - clear it
            for p in already:
                p.unlink(missing_ok=True)

        sra_file = None
        if use_prefetch:
            # Only REUSE a cached .sra that a previous prefetch actually FINISHED, marked by
            # <cache>/.<acc>.sra.done (written only after prefetch exits 0, which includes its
            # own checksum validation). A .sra with no marker is from an interrupted download,
            # so delete it and re-fetch - this is what stops a half-downloaded archive from
            # being split into a truncated FASTQ on resume.
            done = cache_dir / f".{acc}.sra.done"
            cached = _prefetched_sra(cache_dir, acc)
            if cached and done.exists():
                sra_file = cached
            else:
                if cached or done.exists():
                    log.warning("%s: cached .sra is from an interrupted run - removing it "
                                "and re-downloading", acc)
                    _purge_sra_cache(cache_dir, acc)
                    done.unlink(missing_ok=True)
                _clear_stale_locks(cache_dir, acc)   # clear locks left by an interrupted run
                try:
                    get_wrapper("sra", "prefetch")(
                        threads=cfg.threads, log_dir=ctx.logs,
                    ).run(accession=acc, cache_dir=cache_dir,
                          max_size=sra_cfg.get("max_size", "100G"))
                    sra_file = _prefetched_sra(cache_dir, acc)
                    if sra_file:
                        done.touch()          # prefetch exited 0 -> the archive is complete
                except Exception as e:
                    # Not fatal: the dumper can stream straight from NCBI instead.
                    log.warning("prefetch failed for %s (%s) - dumping over the network",
                                acc, e)

        reads_dir.mkdir(parents=True, exist_ok=True)   # create raw_reads/ now dumping begins
        wrapper = get_wrapper("sra", dumper)(threads=cfg.threads, log_dir=ctx.logs)
        kwargs = dict(accession=acc, outdir=reads_dir, split=split, sra_file=sra_file)
        if dumper == "fastq-dump":
            kwargs["gzip"] = want_gzip
        else:
            kwargs["tmp_dir"] = ctx.tmp
        wrapper.run(**kwargs)

        produced = existing_fastqs(reads_dir, acc)
        if not produced:
            raise OutputValidationError(
                f"{dumper} produced no FASTQ for {acc}. Check the accession is a "
                f"run (SRR/ERR/DRR) and that the SRA Toolkit is configured "
                f"(`vdb-config --interactive`)."
            )
        if want_gzip:
            produced = [_gzip_in_place(p) if p.suffix == ".fastq" else p
                        for p in produced]
        _fetch_marker(reads_dir, acc).touch()     # mark this accession complete
        log.info("%s -> %s", acc, ", ".join(p.name for p in produced))

    if sra_cfg.get("keep_sra", False) is False and cache_dir.exists():
        shutil.rmtree(cache_dir, ignore_errors=True)
        cache_dir.mkdir(parents=True, exist_ok=True)

    return reads_dir
