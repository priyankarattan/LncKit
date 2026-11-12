"""Adapter/quality trimming stage (fastp / Trimmomatic / Cutadapt)."""
from __future__ import annotations
from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.trim")


def _complete_gz(path) -> bool:
    """True if `path` is a valid, non-truncated gzip. A trim run killed mid-write leaves
    a truncated .gz; `gzip -t` fails on it, so we re-trim rather than reuse a partial file
    (which would desync R1/R2 pairing downstream)."""
    import subprocess
    try:
        return subprocess.run(["gzip", "-t", str(path)],
                              capture_output=True).returncode == 0
    except Exception:      # gzip missing / unreadable -> treat as incomplete, re-trim
        return False


def run_trim(ctx: RunContext, samples):
    """samples: list of dicts {name, fq1, fq2?}. Returns list of trimmed sample dicts."""
    backend = ctx.config.tools.get("trim", "fastp")
    resume = ctx.config.get("execution.resume", True)
    out = ctx.dirs["trimmed_reads"]
    params = dict(ctx.config.get(f"tool_params.{backend}", {}) or {})
    if backend == "trimmomatic":
        # Trimmomatic is a Java tool with a fixed default heap; give it a heap sized to
        # the run's memory so its parallel gzip compressor doesn't OOM on PE data.
        params.setdefault("mem_gb", ctx.config.get("resources.memory_gb", 8))
    wrapper = get_wrapper("trim", backend)(
        threads=ctx.config.threads,
        extra_params=params,
        log_dir=ctx.logs,
    )
    trimmed = []
    for s in samples:
        o1 = out / f"{s['name']}_1.clean.fastq.gz"
        o2 = out / f"{s['name']}_2.clean.fastq.gz" if s.get("fq2") else None
        # Resume only if the trimmed files EXIST *and* are COMPLETE gzips. A run killed
        # mid-trim (walltime/disconnect) leaves a truncated .gz - reusing it desyncs the
        # R1/R2 pairing (aligners then abort with 'paired reads have different names').
        done = (o1.exists() and _complete_gz(o1)
                and (o2 is None or (o2.exists() and _complete_gz(o2))))
        if done and resume:
            log.info("%s: trimmed reads already present - skipping", s["name"])
        else:
            wrapper.run(
                in1=s["fq1"], out1=o1, in2=s.get("fq2"), out2=o2,
                json_report=ctx.dirs["qc"] / f"{s['name']}.fastp.json",
                html_report=ctx.dirs["qc"] / f"{s['name']}.fastp.html",
            )
        trimmed.append({"name": s["name"], "fq1": o1, "fq2": o2})
    return trimmed
