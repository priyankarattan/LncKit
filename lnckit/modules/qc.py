"""Quality-control stage: FastQC, fastp (report-only), or MultiQC (FastQC + aggregate)."""
from __future__ import annotations
from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.qc")


def run_qc(ctx: RunContext, fastqs):
    backend = ctx.config.tools.get("qc", "fastp")
    outdir = ctx.dirs["qc"]
    threads, logs = ctx.config.threads, ctx.logs

    if backend == "fastqc":
        get_wrapper("qc", "fastqc")(threads=threads, log_dir=logs).run(
            fastqs=list(fastqs), outdir=outdir)

    elif backend == "multiqc":
        # FastQC per file, then MultiQC to roll the reports into one summary.
        get_wrapper("qc", "fastqc")(threads=threads, log_dir=logs).run(
            fastqs=list(fastqs), outdir=outdir)
        get_wrapper("qc", "multiqc")(threads=threads, log_dir=logs).run(
            analysis_dir=outdir, outdir=outdir, name="multiqc_report")

    else:
        # qc=fastp reuses fastp's per-sample HTML/JSON report - but ONLY if fastp is also
        # the trimmer. If a different trimmer is used (trimmomatic/cutadapt), fastp never
        # runs, so nothing would land in qc/. Fall back to FastQC so QC is always produced.
        if ctx.config.tools.get("trim") == "fastp":
            log.info("qc=fastp: the QC report is written by fastp during trimming (qc/)")
        else:
            log.info("qc=fastp but trimmer is %r - running FastQC so a QC report exists",
                     ctx.config.tools.get("trim"))
            get_wrapper("qc", "fastqc")(threads=threads, log_dir=logs).run(
                fastqs=list(fastqs), outdir=outdir)

    return outdir
