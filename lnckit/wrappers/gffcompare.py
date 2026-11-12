"""gffcompare wrapper - compare assembled transcripts to reference annotation.

Produces class codes (u/i/x/o/j/=/...) used to separate novel candidate
transcripts from known/reference transcripts.
"""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("compare", "gffcompare")
class GffcompareWrapper(BaseWrapper):
    requires = ("gffcompare",)
    tool_name = "gffcompare"

    def build_command(self, *, query_gtf: Path, reference_gtf: Path, out_prefix: Path, **_) -> List[str]:
        return ["gffcompare", "-r", str(reference_gtf), "-o", str(out_prefix), str(query_gtf)]

    def expected_outputs(self, *, out_prefix: Path, **_) -> List[Path]:
        return [Path(f"{out_prefix}.tracking")]
