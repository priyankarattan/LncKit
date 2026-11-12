"""CNCI and PLEK wrappers - additional coding-potential backends.

CNCI (Coding-Non-Coding Index) is an alignment-free, model-based classifier; run
with its plant model (`-m pl`) it is one of the two tools (with CPC/CPC2) used in
several plant lncRNA studies. PLEK was trained largely on human/mouse data and
under-performs on plant transcriptomes, so it is offered only for comparison.

Both are Python-2 tools installed from source into the lnckit_py2 env (see
core/installer.py) and exposed on PATH as CNCI.py / PLEK.
"""
from __future__ import annotations
from pathlib import Path
from typing import List

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@register("coding_potential", "cnci")
class CnciWrapper(BaseWrapper):
    requires = ("cnci",)
    tool_name = "cnci"

    def build_command(self, *, fasta: Path, out_dir: Path, model: str = "pl",
                      **_) -> List[str]:
        # -m pl = plant model; CNCI writes CNCI.index into the output directory.
        return ["CNCI.py", "-f", str(fasta), "-o", str(out_dir),
                "-m", model, "-p", self.threads]

    def expected_outputs(self, *, out_dir: Path, **_) -> List[Path]:
        return [Path(out_dir) / "CNCI.index"]


@register("coding_potential", "plek")
class PlekWrapper(BaseWrapper):
    requires = ("plek",)
    tool_name = "plek"
    # PLEK 1.2 exits non-zero even after a fully successful run (a long-standing bug),
    # so success is judged by its output file, not the exit code. See expected_outputs.
    check_returncode = False

    def build_command(self, *, fasta: Path, out_file: Path, **_) -> List[str]:
        # PLEK is exposed on PATH as PLEK.py (the installer's shim name), not "PLEK".
        return ["PLEK.py", "-fasta", str(fasta), "-out", str(out_file),
                "-thread", self.threads]

    def expected_outputs(self, *, out_file: Path, **_) -> List[Path]:
        return [Path(out_file)]
