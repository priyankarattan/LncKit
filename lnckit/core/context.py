"""RunContext: the single object threaded through every step.

Bundles the validated config, the output directory layout, the checkpoint store,
the run metadata, and lazily-resolved tool paths - so individual steps/wrappers
receive one argument and never reach for globals.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

from lnckit.core.checkpoint import CheckpointStore
from lnckit.core.config import Config
from lnckit.core.metadata import RunMetadata

# Canonical output subdirectory layout (mirrors the workflow stages).
SUBDIRS = [
    "logs", "sra_cache", "raw_reads", "merged_reads", "qc", "trimmed_reads",
    "index", "alignment",
    "assembly", "filtered_transcripts", "coding_prediction",
    "protein_filter", "rfam_filter", "orf", "classification", "final",
    "reports", "temporary",
]


class _StepDirs(dict):
    """A name -> Path map that creates each directory the FIRST time it is read, so the
    output tree is built up stage by stage (sra_cache when the download starts, raw_reads
    when dumping starts, alignment when mapping starts, ...) instead of every folder
    appearing at once the moment the run begins."""
    def __getitem__(self, key):
        p = dict.__getitem__(self, key)
        p.mkdir(parents=True, exist_ok=True)
        return p


@dataclass
class RunContext:
    config: Config
    output_dir: Path
    checkpoints: CheckpointStore = field(init=False)
    metadata: RunMetadata = field(init=False)
    dirs: Dict[str, Path] = field(default_factory=_StepDirs)

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # Register each stage's output path but DO NOT create the directory yet - each one
        # is made the first time a stage actually uses it (see _StepDirs), so the output
        # folders appear step by step as the pipeline runs rather than all at once.
        for sub in SUBDIRS:
            dict.__setitem__(self.dirs, sub, self.output_dir / sub)
        self.checkpoints = CheckpointStore(self.output_dir)
        self.metadata = RunMetadata(config=self.config.raw)

    @property
    def tmp(self) -> Path:
        return self.dirs["temporary"]

    @property
    def logs(self) -> Path:
        return self.dirs["logs"]

    @classmethod
    def from_config(cls, config: Config) -> "RunContext":
        return cls(config=config, output_dir=config.output_dir)
