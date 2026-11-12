"""gffread wrapper - extract transcript FASTA from a GTF + genome."""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("extract", "gffread")
class GffreadWrapper(BaseWrapper):
    requires = ("gffread",)
    tool_name = "gffread"

    def build_command(self, *, gtf: Path, genome: Path, fasta_out: Path, **_) -> List[str]:
        return ["gffread", "-w", str(fasta_out), "-g", str(genome), str(gtf)]

    def expected_outputs(self, *, fasta_out: Path, **_) -> List[Path]:
        return [Path(fasta_out)]
