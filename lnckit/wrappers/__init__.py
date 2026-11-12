"""Wrapper classes for external bioinformatics tools.

Design rule: NEVER reimplement an external algorithm. Each wrapper only builds a
command line, executes it through the shared subprocess runner, validates the
declared outputs, and records provenance. New tools are registered in registry.py.

Registration is an import side-effect of the @register decorator, so every wrapper
submodule must actually be imported before the registry is consulted - otherwise
get_wrapper() raises KeyError for tools that plainly exist. Importing them here
means `from lnckit.wrappers import get_wrapper` is always enough.

Third-party plugins register the same way and are discovered through the
"lnckit.wrappers" entry-point group, so a new backend needs no core edit.
"""
import importlib
import pkgutil

from lnckit.wrappers.registry import register, get_wrapper, available  # noqa: F401


def _load_builtin_wrappers() -> None:
    for mod in pkgutil.iter_modules(__path__):
        if mod.name not in {"base", "registry"}:
            importlib.import_module(f"{__name__}.{mod.name}")


def load_plugins(group: str = "lnckit.wrappers") -> list:
    """Import third-party wrapper packages advertised via entry points."""
    loaded = []
    try:
        from importlib.metadata import entry_points
        for ep in entry_points(group=group):
            ep.load()
            loaded.append(ep.name)
    except Exception:  # pragma: no cover - a broken plugin must not kill the run
        pass
    return loaded


_load_builtin_wrappers()
load_plugins()
