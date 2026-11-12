"""Abundance re-estimation against the merged annotation.

Why a second StringTie pass is needed
-------------------------------------
`stringtie --merge` unifies per-sample assemblies into one non-redundant GTF with
new transcript IDs (MSTRG.*), but it carries no expression values - the per-sample
TPMs refer to the old, per-sample IDs. Any TPM filter applied to merged IDs
therefore needs abundances re-estimated *on the merged annotation*, which is what
`stringtie -e -G merged.gtf` does: it quantifies the given transcripts in each BAM
without assembling anything new.

Output: {transcript_id: [TPM in sample 1, TPM in sample 2, ...]}, which lets the
expression filter ask "expressed in at least N samples" rather than collapsing
tissue-specific lncRNAs into a meaningless average.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger
from lnckit.filters.expression import merge_sample_tpms
from lnckit.utils.gtf import transcript_abundance
from lnckit.wrappers import get_wrapper

log = get_logger("lnckit.quantify")


def run_quantify(ctx: RunContext, bams: List[Path], merged_gtf: Path,
                 metric: str = "tpm") -> Dict[str, List[float]]:
    """Re-estimate abundance (TPM or FPKM) for every merged transcript in every sample."""
    if not bams:
        return {}
    outdir = ctx.dirs["assembly"] / "quant"
    outdir.mkdir(parents=True, exist_ok=True)
    strand = ctx.config.get("inputs.strandness")

    per_sample: Dict[str, Dict[str, float]] = {}
    for bam in bams:
        name = Path(bam).stem
        gtf_out = outdir / f"{name}.abund.gtf"
        if not gtf_out.exists() or not ctx.config.get("execution.resume", True):
            get_wrapper("quant", "stringtie")(
                threads=ctx.config.threads, log_dir=ctx.logs,
            ).run(bam=bam, merged_gtf=merged_gtf, gtf_out=gtf_out, strandness=strand)
        per_sample[name] = transcript_abundance(gtf_out, metric)
        log.info("quantified %s: %d transcripts (%s)", name, len(per_sample[name]),
                 metric.upper())

    return merge_sample_tpms(per_sample)
