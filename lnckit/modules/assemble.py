"""Transcript assembly + merge + FASTA extraction (StringTie / Cufflinks).

Per-sample assembly, then a single merged non-redundant annotation across all
samples - a lncRNA expressed in only one tissue must survive into the merged set,
which is why merging happens before, not after, filtering.

`stringtie --merge` is the merge backend for both assemblers: cuffmerge is Python-2
only and effectively unmaintained, and StringTie merges Cufflinks GTFs happily.

The merged GTF plus the genome yields the transcript FASTA (gffread) that every
downstream filter - ORF, coding-potential, protein homology, Rfam - actually reads.
"""
from __future__ import annotations

from lnckit.core.context import RunContext
from lnckit.core.exceptions import ConfigError
from lnckit.core.logging_utils import get_logger
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.assemble")


def run_assemble(ctx: RunContext, bams):
    if not bams:
        log.warning("no BAMs to assemble - skipping assembly")
        return None, None

    backend = ctx.config.tools.get("assembler", "stringtie")
    strand = ctx.config.get("inputs.strandness")
    ann = ctx.config.get("reference.annotation_gtf")
    genome = ctx.config.get("reference.genome_fasta")
    asm_dir = ctx.dirs["assembly"]
    resume = ctx.config.get("execution.resume", True)
    threads, logs = ctx.config.threads, ctx.logs

    per_sample = []
    for bam in bams:
        gtf = asm_dir / f"{bam.stem}.gtf"
        if not (gtf.exists() and resume):
            get_wrapper("assembler", backend)(threads=threads, log_dir=logs).run(
                bam=bam, gtf_out=gtf, annotation=ann, strandness=strand)
        per_sample.append(gtf)

    merged = asm_dir / "merged.gtf"
    if not (merged.exists() and resume):
        get_wrapper("merge", "stringtie")(threads=threads, log_dir=logs).run(
            gtfs=per_sample, merged_out=merged, annotation=ann)

    if not genome:
        raise ConfigError("reference.genome_fasta is required to extract the "
                          "transcript FASTA from the merged assembly.")
    tx_fasta = asm_dir / "transcripts.fasta"
    if not (tx_fasta.exists() and resume):
        get_wrapper("extract", "gffread")(threads=threads, log_dir=logs).run(
            gtf=merged, genome=genome, fasta_out=tx_fasta)

    log.info("assembly: %d per-sample GTF(s) -> %s", len(per_sample), merged)
    return merged, tx_fasta
