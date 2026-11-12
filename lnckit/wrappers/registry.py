"""Plugin registry mapping (stage, backend-name) -> wrapper class.

This is the extension seam: a third-party package can register a new backend
(e.g. a new aligner) with @register("aligner", "minimap2") without editing core
code. The workflow engine only ever asks the registry for a wrapper by name.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Tuple, Type

_REGISTRY: Dict[Tuple[str, str], Type] = {}


def register(stage: str, name: str) -> Callable:
    def deco(cls):
        _REGISTRY[(stage, name)] = cls
        return cls
    return deco


def get_wrapper(stage: str, name: str) -> Type:
    try:
        return _REGISTRY[(stage, name)]
    except KeyError as e:
        raise KeyError(
            f"No wrapper registered for stage='{stage}', backend='{name}'. "
            f"Available: {available(stage)}"
        ) from e


def available(stage: str = None) -> List[str]:
    if stage is None:
        return [f"{s}:{n}" for (s, n) in _REGISTRY]
    return [n for (s, n) in _REGISTRY if s == stage]
