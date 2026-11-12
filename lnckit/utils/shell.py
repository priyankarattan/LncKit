"""Thin, safe wrapper around subprocess used by every tool wrapper.

Centralising subprocess handling means logging, timing, environment control and
error translation happen in exactly one place, so individual wrappers stay tiny.
"""
from __future__ import annotations

import os
import shlex
import subprocess
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from lnckit.core.exceptions import ToolExecutionError
from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.shell")


@dataclass
class RunResult:
    """Outcome of a single external command."""
    cmd: str
    returncode: int
    stdout: str
    stderr: str
    wall_seconds: float
    log_path: Optional[Path] = None
    extra: dict = field(default_factory=dict)


def run(
    argv: Sequence[str],
    tool: str,
    log_file: Optional[Path] = None,
    cwd: Optional[Path] = None,
    env: Optional[dict] = None,
    check: bool = True,
    stdout_path: Optional[Path] = None,
) -> RunResult:
    """Execute *argv*, capturing/streaming output.

    Parameters
    ----------
    argv        : command as a list (never a shell string -> no injection).
    tool        : human-readable tool name for logging/errors.
    log_file    : append combined stderr here (per-step log).
    stdout_path : if given, redirect stdout to this file (for tools that print
                  results to stdout, e.g. diamond, FEELnc_filter).
    check       : raise ToolExecutionError on non-zero exit.
    """
    argv = [str(a) for a in argv]
    printable = " ".join(shlex.quote(a) for a in argv)
    log.info("%s $ %s", tool, printable)   # plain text (markup disabled for robustness)

    full_env = {**os.environ, **(env or {})}
    start = time.time()

    out_handle = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    try:
        proc = subprocess.run(
            argv,
            cwd=str(cwd) if cwd else None,
            env=full_env,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    finally:
        if stdout_path:
            out_handle.close()

    wall = time.time() - start
    stdout = "" if stdout_path else (proc.stdout or "")
    stderr = proc.stderr or ""

    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as fh:
            fh.write(f"\n$ {printable}\n{stderr}\n[exit {proc.returncode} | {wall:.1f}s]\n")

    if check and proc.returncode != 0:
        raise ToolExecutionError(tool, proc.returncode, printable, stderr)

    log.debug("%s finished in %.1fs (exit %d)", tool, wall, proc.returncode)
    return RunResult(printable, proc.returncode, stdout, stderr, wall, log_file)


def run_pipeline(
    producer: Sequence[str],
    consumer: Sequence[str],
    tool: str,
    producer_log: Optional[Path] = None,
    consumer_log: Optional[Path] = None,
    env: Optional[dict] = None,
    check: bool = True,
) -> RunResult:
    """Run `producer | consumer`, streaming the producer's stdout into the consumer's stdin.

    Used to pipe an aligner (which emits SAM on stdout) straight into `samtools sort`,
    so the multi-GB intermediate SAM never touches disk - the single biggest cause of
    the alignment stage filling the filesystem on a large library. Both processes run
    concurrently; each stderr is drained in its own thread so a full stderr pipe can
    never deadlock the stream.
    """
    prod = [str(a) for a in producer]
    cons = [str(a) for a in consumer]
    prod_str = " ".join(shlex.quote(a) for a in prod)
    cons_str = " ".join(shlex.quote(a) for a in cons)
    printable = f"{prod_str} | {cons_str}"
    log.info("%s $ %s", tool, printable)

    full_env = {**os.environ, **(env or {})}
    start = time.time()

    p1 = subprocess.Popen(prod, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          env=full_env, text=True)
    p2 = subprocess.Popen(cons, stdin=p1.stdout, stdout=subprocess.DEVNULL,
                          stderr=subprocess.PIPE, env=full_env, text=True)
    # Close our copy of the read end so the producer gets SIGPIPE if the consumer exits.
    p1.stdout.close()

    collected: dict = {}

    def _drain(proc, key):
        collected[key] = proc.stderr.read()
        proc.stderr.close()

    t1 = threading.Thread(target=_drain, args=(p1, "prod"))
    t2 = threading.Thread(target=_drain, args=(p2, "cons"))
    t1.start()
    t2.start()
    p2.wait()
    p1.wait()
    t1.join()
    t2.join()
    wall = time.time() - start

    prod_err = collected.get("prod") or ""
    cons_err = collected.get("cons") or ""

    if producer_log:
        Path(producer_log).parent.mkdir(parents=True, exist_ok=True)
        with open(producer_log, "a") as fh:
            fh.write(f"\n$ {prod_str}\n{prod_err}\n[exit {p1.returncode}]\n")
    if consumer_log:
        Path(consumer_log).parent.mkdir(parents=True, exist_ok=True)
        with open(consumer_log, "a") as fh:
            fh.write(f"\n$ {cons_str}\n{cons_err}\n[exit {p2.returncode} | {wall:.1f}s]\n")

    if check and p1.returncode != 0:
        raise ToolExecutionError(tool, p1.returncode, prod_str, prod_err)
    if check and p2.returncode != 0:
        raise ToolExecutionError(tool, p2.returncode, cons_str, cons_err)

    log.debug("%s pipeline finished in %.1fs (producer exit %s, consumer exit %s)",
              tool, wall, p1.returncode, p2.returncode)
    return RunResult(printable, p2.returncode, "", prod_err + cons_err, wall, consumer_log)
