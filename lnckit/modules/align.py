"""Alignment + BAM sorting stage (HISAT2 / STAR / Bowtie2 / BWA + samtools).

STAR writes a coordinate-sorted BAM itself. Every other backend emits SAM: by
default LncKit streams that SAM straight into `samtools sort` through a pipe, so
the coordinate-sorted BAM is the only thing written to disk. The uncompressed SAM
for one RNA-seq library is several times the size of the BAM (tens of GB) and
nothing downstream reads it, so never materialising it is what keeps a large run
from filling the filesystem mid-alignment. Set execution.keep_intermediate to
write the SAM to disk anyway (then sort it); a SAM already on disk is reused as-is.

Only splice-aware aligners (HISAT2/STAR) are valid ahead of transcript assembly;
core.config.validate() rejects the other combinations before any read is mapped.
"""
from __future__ import annotations

from pathlib import Path

from lnckit.core.context import RunContext
from lnckit.core.exceptions import OutputValidationError
from lnckit.core.logging_utils import get_logger
from lnckit.utils.shell import run_pipeline
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.align")


def _stream_align_to_bam(ctx, backend, index_ref, bam, sample,
                         strand, dta, threads, logs):
    """Align and coordinate-sort in a single pass.

    The aligner writes SAM to stdout, piped straight into `samtools sort`, so the
    multi-GB intermediate SAM is never written to disk - this is what stops the
    alignment stage from filling the filesystem on a large library. The BAM is
    written to a `.part.bam` first and renamed only on success, so an interrupted
    run never leaves a truncated BAM that resume would mistake for a finished one.
    """
    aligner = get_wrapper("aligner", backend)(threads=threads, log_dir=logs)
    sorter = get_wrapper("bam", "samtools")(threads=threads, log_dir=logs)
    aligner.check_dependencies()
    sorter.check_dependencies()

    bam_tmp = bam.with_suffix(".part.bam")
    bam_tmp.unlink(missing_ok=True)

    producer = aligner.build_command(
        index_prefix=index_ref, sam_out=None,
        fq1=sample["fq1"], fq2=sample.get("fq2"), strandness=strand, dta=dta)
    consumer = sorter.build_command(sam_in=Path("-"), bam_out=bam_tmp)
    run_pipeline(producer, consumer, tool="%s | samtools-sort" % backend,
                 producer_log=logs / ("%s.log" % aligner.tool_name),
                 consumer_log=logs / ("%s.log" % sorter.tool_name),
                 env=aligner.env())

    if not bam_tmp.exists() or bam_tmp.stat().st_size == 0:
        bam_tmp.unlink(missing_ok=True)
        raise OutputValidationError(
            "%s | samtools-sort: no sorted BAM produced for %s" % (backend, bam.name))
    bam_tmp.replace(bam)


def _index_bam(ctx: RunContext, bam, threads, logs):
    """Create <bam>.bai for a coordinate-sorted BAM (idempotent on resume)."""
    if (bam.with_suffix(bam.suffix + ".bai").exists()
            and ctx.config.get("execution.resume", True)):
        return
    get_wrapper("bam_index", "samtools")(threads=threads, log_dir=logs).run(bam_in=bam)


def _discover_alignment_inputs(aln_dir):
    """Build a name-only sample list from SAM/BAM already in the alignment dir.

    Lets a run JUMP straight to alignment (lnckit.sh --from align) when you already have
    SAM/BAM files but no reads - each `<name>.sam`/`<name>.bam` becomes one sample.
    """
    names = {p.name[:-4] for p in aln_dir.glob("*.sam")}
    names |= {p.stem for p in aln_dir.glob("*.bam")}
    return [{"name": n, "fq1": None, "fq2": None} for n in sorted(names)]


def run_align(ctx: RunContext, samples, index_ref):
    backend = ctx.config.tools.get("aligner", "hisat2")
    strand = ctx.config.get("inputs.strandness")
    if strand == "auto":                 # still unresolved -> align unstranded; the engine
        strand = None                    # detects RF/FR from these BAMs before assembly
    # HISAT2 --dta: on by default (recommended for StringTie). Set tool_params.hisat2.dta
    # to false to match pipelines that omit it (e.g. Plant-LncPipe, for exact replication).
    dta = ctx.config.get("tool_params.hisat2.dta", True)
    aln_dir = ctx.dirs["alignment"]
    resume = ctx.config.get("execution.resume", True)
    keep_sam = ctx.config.get("execution.keep_intermediate", False)
    threads, logs = ctx.config.threads, ctx.logs
    bams = []

    # Jumped straight to align (no upstream samples)? Use SAM/BAM already on disk.
    if not samples:
        samples = _discover_alignment_inputs(aln_dir)
        if samples:
            log.info("no upstream samples - found %d alignment file(s) to process",
                     len(samples))

    for s in samples:
        name = s["name"]
        if backend == "star":
            bam = aln_dir / f"{name}.Aligned.sortedByCoord.out.bam"
            if not (bam.exists() and resume):
                get_wrapper("aligner", "star")(threads=threads, log_dir=logs).run(
                    index_dir=index_ref, out_prefix=aln_dir / f"{name}.",
                    fq1=s["fq1"], fq2=s.get("fq2"))
            _index_bam(ctx, bam, threads, logs)         # STAR BAM is sorted -> index it
            bams.append(bam)
            continue

        bam = aln_dir / f"{name}.bam"
        if bam.exists() and resume:
            log.info("%s: sorted BAM already present - skipping alignment", name)
            _index_bam(ctx, bam, threads, logs)         # ensure the .bai exists
            bams.append(bam)
            continue

        # Default path: stream the aligner straight into samtools sort so no SAM is
        # ever written to disk (the SAM for one RNA-seq library is tens of GB and is
        # what fills the filesystem). Two exceptions keep the old behaviour:
        #   * a SAM is already on disk (you supplied it, or a prior run stopped after
        #     alignment) -> reuse it and only do the SAM->BAM conversion;
        #   * keep_intermediate is set -> write the SAM, then sort it.
        sam = aln_dir / f"{name}.sam"
        if sam.exists() and resume:
            log.info("%s: reusing existing SAM %s - skipping the aligner", name, sam.name)
            get_wrapper("bam", "samtools")(threads=threads, log_dir=logs).run(
                sam_in=sam, bam_out=bam)
            if not keep_sam:
                sam.unlink(missing_ok=True)
        elif s.get("fq1"):
            if keep_sam:
                get_wrapper("aligner", backend)(threads=threads, log_dir=logs).run(
                    index_prefix=index_ref, sam_out=sam,
                    fq1=s["fq1"], fq2=s.get("fq2"), strandness=strand, dta=dta)
                get_wrapper("bam", "samtools")(threads=threads, log_dir=logs).run(
                    sam_in=sam, bam_out=bam)
            else:
                _stream_align_to_bam(ctx, backend, index_ref, bam, s,
                                     strand, dta, threads, logs)
        else:
            log.warning("%s: no reads and no SAM/BAM - cannot align; skipping", name)
            continue
        _index_bam(ctx, bam, threads, logs)
        log.info("%s: sorted BAM %s (+ .bai index)", name, bam.name)
        bams.append(bam)

    return bams
