"""Transcript-to-reference comparison stage (gffcompare) + novel extraction."""
from __future__ import annotations
from lnckit.core.context import RunContext
from lnckit.wrappers import get_wrapper


def run_compare(ctx: RunContext, merged_gtf):
    ann = ctx.config.get("reference.annotation_gtf")
    # gffcompare output is a pure intermediate (feeds the identify step), so it lives
    # under temporary/ rather than as its own results folder.
    cmp_dir = ctx.tmp / "gffcompare"
    cmp_dir.mkdir(parents=True, exist_ok=True)
    prefix = cmp_dir / "gffcmp"
    get_wrapper("compare", "gffcompare")(threads=ctx.config.threads, log_dir=ctx.logs).run(
        query_gtf=merged_gtf, reference_gtf=ann, out_prefix=prefix)
    return prefix
