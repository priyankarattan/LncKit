"""samtools wrapper - SAM->sorted BAM conversion and indexing."""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("bam", "samtools")
class SamtoolsSortWrapper(BaseWrapper):
    requires = ("samtools",)
    tool_name = "samtools-sort"

    def build_command(self, *, sam_in: Path, bam_out: Path, **_) -> List[str]:
        # samtools sort accepts SAM/BAM directly; sorts into a coordinate-sorted BAM.
        return ["samtools", "sort", "-@", self.threads, "-o", str(bam_out), str(sam_in)]

    def expected_outputs(self, *, bam_out: Path, **_) -> List[Path]:
        return [Path(bam_out)]


@register("bam_index", "samtools")
class SamtoolsIndexWrapper(BaseWrapper):
    requires = ("samtools",)
    tool_name = "samtools-index"

    def build_command(self, *, bam_in: Path, **_) -> List[str]:
        # Index a coordinate-sorted BAM -> <bam>.bai (needed by StringTie -e, IGV, etc.).
        return ["samtools", "index", "-@", self.threads, str(bam_in)]

    def expected_outputs(self, *, bam_in: Path, **_) -> List[Path]:
        return [Path(f"{bam_in}.bai")]
