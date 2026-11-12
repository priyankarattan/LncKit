"""Run-level metadata / provenance record.

Reproducibility is a first-class goal: every run writes run_metadata.json capturing
the exact config, resolved tool paths+versions, per-step runtime, and transcript
counts after each filter. This file backs the HTML/PDF report and lets a reviewer
reconstruct precisely what was executed.
"""
from __future__ import annotations

import json
import platform
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List

from lnckit.__about__ import __version__


def _platform_string() -> str:
    return platform.platform()


def _python_version() -> str:
    return platform.python_version()


@dataclass
class StepRecord:
    name: str
    tool: str = ""
    status: str = "pending"
    runtime_seconds: float = 0.0
    transcripts_in: int = -1
    transcripts_out: int = -1
    outputs: List[str] = field(default_factory=list)


@dataclass
class RunMetadata:
    lnckit_version: str = __version__
    started: float = field(default_factory=time.time)
    finished: float = 0.0
    platform: str = field(default_factory=_platform_string)
    python: str = field(default_factory=_python_version)
    config: dict = field(default_factory=dict)
    tool_versions: Dict[str, str] = field(default_factory=dict)
    steps: List[StepRecord] = field(default_factory=list)
    final_counts: Dict[str, int] = field(default_factory=dict)

    def add_step(self, rec: StepRecord) -> None:
        self.steps.append(rec)

    def write(self, output_dir: Path) -> Path:
        self.finished = time.time()
        out = Path(output_dir) / "reports" / "run_metadata.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        out.write_text(json.dumps(payload, indent=2, default=str))
        return out
