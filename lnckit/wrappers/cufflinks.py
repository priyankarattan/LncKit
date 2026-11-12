"""Cufflinks wrapper - alternative transcript assembler to StringTie.

Cufflinks writes a fixed filename (`transcripts.gtf`) into an output *directory*,
so the wrapper runs it in a per-sample directory and then moves the result to the
`gtf_out` path the assembly module expects - keeping the assembler backends
interchangeable behind one interface.

Note: `cuffmerge` is Python-2 only and effectively unmaintained, so merging is
delegated to `stringtie --merge` (which happily merges Cufflinks GTFs). That is
also what current Cufflinks-based pipelines do in practice.
"""
from __future__ import annotations
import shutil
from pathlib import Path
from typing import List, Optional

from lnckit.utils.shell import RunResult
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register

_LIBTYPE = {"RF": "fr-firststrand", "FR": "fr-secondstrand"}


@register("assembler", "cufflinks")
class CufflinksWrapper(BaseWrapper):
    requires = ("cufflinks",)
    tool_name = "cufflinks"

    def _workdir(self, gtf_out: Path) -> Path:
        return Path(gtf_out).parent / f"{Path(gtf_out).stem}_cufflinks"

    def build_command(self, *, bam: Path, gtf_out: Path,
                      annotation: Optional[Path] = None,
                      strandness: Optional[str] = None, **_) -> List[str]:
        workdir = self._workdir(gtf_out)
        workdir.mkdir(parents=True, exist_ok=True)
        cmd = ["cufflinks", "-p", self.threads, "-o", str(workdir)]
        if strandness in _LIBTYPE:
            cmd += ["--library-type", _LIBTYPE[strandness]]
        if annotation:
            # -g = RABT guided assembly (keeps novel transcripts, unlike -G)
            cmd += ["-g", str(annotation)]
        cmd.append(str(bam))
        return cmd

    def expected_outputs(self, *, gtf_out: Path, **_) -> List[Path]:
        return [self._workdir(gtf_out) / "transcripts.gtf"]

    def run(self, **kwargs) -> RunResult:
        result = super().run(**kwargs)
        gtf_out = Path(kwargs["gtf_out"])
        shutil.move(str(self._workdir(gtf_out) / "transcripts.gtf"), str(gtf_out))
        return result
