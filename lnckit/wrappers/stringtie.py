"""StringTie wrapper - reference-guided transcript assembly and merge."""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("assembler", "stringtie")
class StringtieAssembleWrapper(BaseWrapper):
    requires = ("stringtie",)
    tool_name = "stringtie"

    def build_command(self, *, bam: Path, gtf_out: Path,
                      annotation: Optional[Path] = None,
                      strandness: Optional[str] = None, **_) -> List[str]:
        cmd = ["stringtie", "-p", self.threads]
        if strandness == "RF":
            cmd += ["--rf"]
        elif strandness == "FR":
            cmd += ["--fr"]
        if annotation:
            cmd += ["-G", str(annotation)]
        cmd += ["-o", str(gtf_out), str(bam)]
        return cmd

    def expected_outputs(self, *, gtf_out: Path, **_) -> List[Path]:
        return [Path(gtf_out)]


@register("quant", "stringtie")
class StringtieQuantWrapper(BaseWrapper):
    """Re-estimate abundance against the *merged* annotation (`stringtie -e`).

    The merged GTF carries no TPM: StringTie's per-sample assemblies use their own
    transcript IDs, and `--merge` does not propagate expression. A TPM filter on the
    merged IDs therefore needs this second, cheap pass, which quantifies every
    merged transcript (MSTRG.*) in every sample. `-e` restricts StringTie to the
    given annotation, so no new transcripts are invented here.
    """

    requires = ("stringtie",)
    tool_name = "stringtie-quant"

    def build_command(self, *, bam: Path, merged_gtf: Path, gtf_out: Path,
                      strandness: Optional[str] = None, **_) -> List[str]:
        cmd = ["stringtie", "-e", "-B", "-p", self.threads, "-G", str(merged_gtf)]
        if strandness == "RF":
            cmd += ["--rf"]
        elif strandness == "FR":
            cmd += ["--fr"]
        cmd += ["-o", str(gtf_out), str(bam)]
        return cmd

    def expected_outputs(self, *, gtf_out: Path, **_) -> List[Path]:
        return [Path(gtf_out)]


@register("merge", "stringtie")
class StringtieMergeWrapper(BaseWrapper):
    requires = ("stringtie",)
    tool_name = "stringtie-merge"

    def build_command(self, *, gtfs: List[Path], merged_out: Path,
                      annotation: Optional[Path] = None, **_) -> List[str]:
        cmd = ["stringtie", "--merge", "-p", self.threads, "-o", str(merged_out)]
        if annotation:
            cmd += ["-G", str(annotation)]
        cmd += [str(g) for g in gtfs]
        return cmd

    def expected_outputs(self, *, merged_out: Path, **_) -> List[Path]:
        return [Path(merged_out)]
