"""Abstract base class shared by every external-tool wrapper.

Concrete wrappers implement two small methods:
  * build_command() -> argv list
  * expected_outputs() -> list[Path]
The base class provides dependency checking, execution, timing, output
validation, checkpoint creation and metadata recording - so each wrapper stays
a few lines and behaves identically across the pipeline.
"""
from __future__ import annotations

import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from lnckit.core.checkpoint import Checkpoint, CheckpointStore
from lnckit.core.dependencies import ResolvedTool, resolve_tool
from lnckit.core.exceptions import OutputValidationError
from lnckit.core.logging_utils import get_logger
from lnckit.utils.shell import RunResult, run


class BaseWrapper(ABC):
    """Base for all tool wrappers."""

    #: canonical dependency keys required by this wrapper (see dependencies.KNOWN_TOOLS)
    requires: Sequence[str] = ()
    #: the tool name reported in logs/metadata
    tool_name: str = "tool"
    #: Search tools (DIAMOND/BLAST/cmscan) legitimately produce a zero-hit table:
    #: a lncRNA candidate set with no protein homology is the *expected* outcome,
    #: not a failure, so those wrappers opt out of the non-empty output check.
    allow_empty_outputs: bool = False
    #: Some tools (notably PLEK) exit non-zero even on a fully successful run. Such
    #: wrappers set this False and are judged solely by expected_outputs() instead -
    #: a genuine failure still surfaces there as a missing/empty output file.
    check_returncode: bool = True

    def __init__(self, threads: int = 4, extra_params: Optional[dict] = None,
                 log_dir: Optional[Path] = None):
        self.threads = threads
        self.extra_params = extra_params or {}
        self.log_dir = Path(log_dir) if log_dir else None
        self.log = get_logger(f"lnckit.wrappers.{self.tool_name}")
        self._resolved: Dict[str, ResolvedTool] = {}

    # ---- to be implemented by subclasses ---------------------------------
    @abstractmethod
    def build_command(self, **kwargs) -> List[str]:
        """Return the argv list for the external tool."""

    @abstractmethod
    def expected_outputs(self, **kwargs) -> List[Path]:
        """Return the files this invocation must produce."""

    # ---- shared machinery ------------------------------------------------
    def check_dependencies(self) -> Dict[str, ResolvedTool]:
        for key in self.requires:
            self._resolved[key] = resolve_tool(key)
        return self._resolved

    def versions(self) -> Dict[str, str]:
        return {k: (rt.version or "unknown") for k, rt in self._resolved.items()}

    def _validate_outputs(self, outputs: List[Path]) -> None:
        for out in outputs:
            p = Path(out)
            if not p.exists():
                raise OutputValidationError(f"{self.tool_name}: expected output missing: {p}")
            if p.is_file() and p.stat().st_size == 0 and not self.allow_empty_outputs:
                raise OutputValidationError(f"{self.tool_name}: output is empty: {p}")

    def env(self, **kwargs) -> Optional[dict]:
        """Extra environment variables for the subprocess (merged over os.environ).

        Default None = inherit the parent environment unchanged. Wrappers whose tool
        needs a specific interpreter/library environment override this - e.g. FEELnc,
        which must run under its own isolated Perl rather than whatever is first on PATH.
        """
        return None

    def run(self, stdout_path: Optional[Path] = None, **kwargs) -> RunResult:
        """Full lifecycle: deps -> command -> execute -> validate outputs."""
        self.check_dependencies()
        argv = self.build_command(**kwargs)
        log_file = (self.log_dir / f"{self.tool_name}.log") if self.log_dir else None
        start = time.time()
        result = run(argv, tool=self.tool_name, log_file=log_file,
                     stdout_path=stdout_path, env=self.env(**kwargs),
                     check=self.check_returncode)
        result.wall_seconds = round(time.time() - start, 2)
        self._validate_outputs(self.expected_outputs(**kwargs))
        result.extra["versions"] = self.versions()
        return result

    def make_checkpoint(self, step: str, started: float, inputs, outputs) -> Checkpoint:
        return CheckpointStore.make(
            step=step, started=started, inputs=inputs, outputs=outputs,
            versions=self.versions(),
        )
