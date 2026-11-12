"""fastp wrapper - QC + adapter/quality trimming in one pass."""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("qc", "fastp")
@register("trim", "fastp")
class FastpWrapper(BaseWrapper):
    requires = ("fastp",)
    tool_name = "fastp"

    def build_command(self, *, in1: Path, out1: Path,
                      in2: Optional[Path] = None, out2: Optional[Path] = None,
                      json_report: Optional[Path] = None,
                      html_report: Optional[Path] = None, **_) -> List[str]:
        cmd = ["fastp", "-i", in1, "-o", out1, "-w", min(self.threads, 16)]
        if in2 and out2:
            cmd += ["-I", in2, "-O", out2]
        if json_report:
            cmd += ["-j", json_report]
        if html_report:
            cmd += ["-h", html_report]
        for k, v in self.extra_params.items():
            cmd += [f"--{k}"] + ([str(v)] if v not in (True, None) else [])
        return cmd

    def expected_outputs(self, *, out1: Path, out2: Optional[Path] = None, **_) -> List[Path]:
        return [out1] + ([out2] if out2 else [])
