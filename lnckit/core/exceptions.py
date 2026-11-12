"""Exception hierarchy for lnckit.

A single well-defined hierarchy lets the CLI translate any failure into a clean
message + exit code, and lets the workflow engine decide what is recoverable
(e.g. resume) versus fatal (e.g. a missing dependency).
"""
from __future__ import annotations


class LncKitError(Exception):
    """Base class for all lnckit errors."""


class ConfigError(LncKitError):
    """Invalid or inconsistent configuration."""


class DependencyError(LncKitError):
    """A required external executable is missing or the wrong version."""


class OutdatedToolError(DependencyError):
    """The executable is present but too old to be used.

    Distinguished from a plain DependencyError so the doctor can say "upgrade this"
    instead of "install this" - telling a user to install a tool they already have
    is how they end up with two copies and a confusing PATH.
    """

    def __init__(self, message: str, tool: str = "", path: str = "", version: str = ""):
        self.tool = tool
        self.path = path
        self.version = version
        super().__init__(message)


class BrokenToolError(DependencyError):
    """The executable exists on PATH but cannot run (exit 126/127).

    Typically a missing shared library or a wrong architecture. This is worse than
    a missing tool, because a naive `shutil.which()` check passes and the pipeline
    then dies in the middle of a long run - so it is detected up front.
    """

    def __init__(self, message: str, tool: str = "", path: str = ""):
        self.tool = tool
        self.path = path
        super().__init__(message)


class ToolExecutionError(LncKitError):
    """An external tool exited with a non-zero status."""

    def __init__(self, tool: str, returncode: int, cmd: str, stderr: str = ""):
        self.tool = tool
        self.returncode = returncode
        self.cmd = cmd
        self.stderr = stderr
        super().__init__(
            f"{tool} failed (exit {returncode}).\n  command: {cmd}\n  stderr: {stderr[-2000:]}"
        )


class OutputValidationError(LncKitError):
    """A tool completed but its declared outputs are missing/empty/malformed."""


class CompatibilityError(ConfigError):
    """Two selected tools are incompatible (e.g. non-splice-aware aligner + StringTie)."""


class CheckpointError(LncKitError):
    """Checkpoint state is corrupt or inconsistent."""
