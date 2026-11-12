"""CPAT wrapper - alignment-free coding-potential assessment.

Uses the plant-retrained hexamer table + logit model (CPAT-plant). The plant
coding-probability cutoff (0.46) is applied downstream in filters/coding_potential.
"""
from __future__ import annotations
import shutil
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


def _cpat_exe() -> str:
    """Return the CPAT executable name that's actually on PATH.

    CPAT's command name changed across versions: the old (Python-2) CPAT and some v3
    builds ship 'cpat.py', while newer v3 (>=3.0.x) installs it as plain 'cpat'. Prefer
    'cpat.py' for back-compat, fall back to 'cpat'; default to 'cpat.py' so the error
    message still names something sensible if neither is found.
    """
    return shutil.which("cpat.py") or shutil.which("cpat") or "cpat.py"


@register("coding_potential", "cpat")
class CpatWrapper(BaseWrapper):
    requires = ("cpat",)
    tool_name = "cpat"

    def build_command(self, *, fasta: Path, out_prefix: Path,
                      hexamer: Path, logit_model: Path,
                      executable: str = None, **_) -> List[str]:
        # `executable` overrides the auto-resolved command - e.g. a CPAT 1.2.4 (Python-2)
        # shim for exact Plant-LncPipe replication. Otherwise use cpat.py/cpat on PATH (v3).
        # Both versions share the -x/-d/-g/-o CLI, and parse_cpat reads either's output.
        return [executable or _cpat_exe(), "-x", str(hexamer), "-d", str(logit_model),
                "-g", str(fasta), "-o", str(out_prefix)]

    def expected_outputs(self, *, out_prefix: Path, **_) -> List[Path]:
        # CPAT writes <prefix> or <prefix>.ORF_prob.tsv depending on version; the
        # coding-potential filter tolerates either.
        return [Path(str(out_prefix)).parent]
