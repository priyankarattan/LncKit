"""Filter protocol + a small result container."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Set


@dataclass
class FilterResult:
    name: str
    kept: Set[str]
    removed: Set[str] = field(default_factory=set)

    @property
    def n_in(self) -> int:
        return len(self.kept) + len(self.removed)

    @property
    def n_out(self) -> int:
        return len(self.kept)
