"""CPC2 wrapper - species-neutral coding-potential classifier."""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("coding_potential", "cpc2")
class Cpc2Wrapper(BaseWrapper):
    requires = ("cpc2",)
    tool_name = "cpc2"

    def build_command(self, *, fasta: Path, out_prefix: Path, **_) -> List[str]:
        return ["CPC2.py", "-i", str(fasta), "-o", str(out_prefix)]

    def expected_outputs(self, *, out_prefix: Path, **_) -> List[Path]:
        return [Path(f"{out_prefix}.txt")]
