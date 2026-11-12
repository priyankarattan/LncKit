"""FastQC wrapper - read-level quality report (no trimming)."""
from __future__ import annotations
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("qc", "fastqc")
class FastQCWrapper(BaseWrapper):
    requires = ("fastqc",)
    tool_name = "fastqc"

    def build_command(self, *, fastqs: List[Path], outdir: Path, **_) -> List[str]:
        return ["fastqc", "-t", self.threads, "-o", outdir, *[str(f) for f in fastqs]]

    def expected_outputs(self, *, fastqs: List[Path], outdir: Path, **_) -> List[Path]:
        # Validate the ACTUAL per-file reports FastQC writes (<name>_fastqc.zip), not just
        # the output directory - otherwise a silent FastQC failure passes as "done" and the
        # qc/ folder ends up empty.
        out = Path(outdir)
        reports = []
        for f in fastqs:
            base = Path(f).name
            for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                if base.endswith(suf):
                    base = base[: -len(suf)]
                    break
            reports.append(out / f"{base}_fastqc.zip")
        return reports or [out]
