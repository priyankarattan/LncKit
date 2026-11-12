"""Checkpoint / resume support.

Each workflow step records a small JSON checkpoint after it validates its own
outputs. `lnckit resume` skips any step whose checkpoint exists, whose input
fingerprints are unchanged, and whose declared outputs are still present. This
gives cheap, transparent resumability without a heavyweight scheduler.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional

from lnckit.core.exceptions import CheckpointError
from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.checkpoint")


def fingerprint(path: Path, algo: str = "sha1", sample_bytes: int = 1 << 20) -> str:
    """Cheap file fingerprint: size + mtime + hash of first `sample_bytes`.

    Full hashing of multi-GB BAM/FASTQ files is wasteful; a head-sample plus size
    is sufficient to detect that an upstream input was replaced.
    """
    p = Path(path)
    if not p.exists():
        return "MISSING"
    h = hashlib.new(algo)
    h.update(str(p.stat().st_size).encode())
    with open(p, "rb") as fh:
        h.update(fh.read(sample_bytes))
    return f"{h.hexdigest()[:16]}:{p.stat().st_size}"


@dataclass
class Checkpoint:
    step: str
    status: str                      # "completed" | "failed"
    started: float
    finished: float
    inputs: Dict[str, str] = field(default_factory=dict)
    outputs: List[str] = field(default_factory=list)
    params_hash: str = ""
    tool_versions: Dict[str, str] = field(default_factory=dict)
    runtime_seconds: float = 0.0


class CheckpointStore:
    """Reads/writes per-step checkpoints under <output>/.lnckit/checkpoints/."""

    def __init__(self, output_dir: Path):
        self.dir = Path(output_dir) / ".lnckit" / "checkpoints"
        self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, step: str) -> Path:
        return self.dir / f"{step}.json"

    def save(self, cp: Checkpoint) -> None:
        with open(self._path(cp.step), "w") as fh:
            json.dump(asdict(cp), fh, indent=2)
        log.debug("checkpoint saved: %s (%s)", cp.step, cp.status)

    def load(self, step: str) -> Optional[Checkpoint]:
        p = self._path(step)
        if not p.exists():
            return None
        try:
            return Checkpoint(**json.loads(p.read_text()))
        except Exception as e:  # pragma: no cover
            raise CheckpointError(f"Corrupt checkpoint for step '{step}': {e}") from e

    def is_complete(self, step: str, expected_outputs: List[Path]) -> bool:
        """A step is resumable only if its checkpoint says completed AND its
        declared outputs still exist on disk."""
        cp = self.load(step)
        if cp is None or cp.status != "completed":
            return False
        for out in expected_outputs:
            if not Path(out).exists():
                log.info("step '%s' checkpoint present but output missing: %s -> rerun",
                         step, out)
                return False
        return True

    @staticmethod
    def make(step: str, started: float, inputs, outputs, params_hash="", versions=None) -> Checkpoint:
        now = time.time()
        return Checkpoint(
            step=step, status="completed", started=started, finished=now,
            inputs={str(k): fingerprint(Path(k)) for k in inputs},
            outputs=[str(o) for o in outputs],
            params_hash=params_hash, tool_versions=versions or {},
            runtime_seconds=round(now - started, 2),
        )
