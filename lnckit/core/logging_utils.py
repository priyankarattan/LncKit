"""Centralised logging.

One configuration point so that every module, wrapper and CLI command logs the
same way: a coloured console stream (via rich) plus a plain-text run log written
into <output>/logs/lnckit.log for reproducibility and post-mortem debugging.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

try:
    from rich.logging import RichHandler
    _HAVE_RICH = True
except Exception:  # pragma: no cover
    _HAVE_RICH = False

_CONFIGURED = False


def setup_logging(
    log_dir: Optional[Path] = None,
    level: str = "INFO",
    quiet: bool = False,
) -> logging.Logger:
    """Configure the root 'lnckit' logger. Idempotent."""
    global _CONFIGURED
    logger = logging.getLogger("lnckit")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    if not quiet:
        if _HAVE_RICH:
            # markup=False: log messages routinely contain file paths and error text with
            # '[...]' (e.g. a '[/.../SRR.sra.lock]' path), which rich would try to parse as
            # markup tags and crash on (MarkupError) - masking the real error. Plain text
            # is safe; the tiny cosmetic markup we used is dropped in favour of robustness.
            console = RichHandler(rich_tracebacks=True, show_path=False, markup=False)
            console.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
        else:  # pragma: no cover
            console = logging.StreamHandler()
            console.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(console)

    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / "lnckit.log")
        fh.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)-8s | %(name)s | %(message)s")
        )
        logger.addHandler(fh)

    logger.propagate = False
    _CONFIGURED = True
    return logger


def get_logger(name: str = "lnckit") -> logging.Logger:
    """Return a namespaced child logger (e.g. get_logger('lnckit.wrappers.hisat2'))."""
    if not _CONFIGURED:
        setup_logging()
    return logging.getLogger(name)
