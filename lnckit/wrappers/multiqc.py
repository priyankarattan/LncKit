"""MultiQC wrapper - aggregate many per-sample QC reports into one summary.

MultiQC is not a per-sample QC tool like FastQC; it *scans a directory* for reports
other tools already wrote (FastQC zips, fastp JSON, HISAT2/STAR logs, ...) and rolls
them into a single interactive HTML. So the `multiqc` QC backend runs FastQC first to
produce the per-sample reports, then MultiQC to combine them.
"""
from __future__ import annotations
from pathlib import Path
from typing import List

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("qc", "multiqc")
class MultiQCWrapper(BaseWrapper):
    requires = ("multiqc",)
    tool_name = "multiqc"

    def build_command(self, *, analysis_dir: Path, outdir: Path,
                      name: str = "multiqc_report", **_) -> List[str]:
        return ["multiqc", str(analysis_dir), "-o", str(outdir),
                "-n", name, "-f"]        # -f overwrites a previous report

    def expected_outputs(self, *, outdir: Path, name: str = "multiqc_report",
                         **_) -> List[Path]:
        return [Path(outdir) / f"{name}.html"]
