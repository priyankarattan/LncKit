"""External-tool discovery and version checking.

lnckit integrates existing bioinformatics software rather than reimplementing
it, so the single most common failure mode is a missing/incompatible executable.
This module resolves each tool once, records its path + version into run
metadata, and fails fast with an actionable message before any long step starts.
"""
from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from typing import Optional

from packaging import version as _v

from lnckit.core.exceptions import (
    BrokenToolError, DependencyError, OutdatedToolError)
from lnckit.core.logging_utils import get_logger
from lnckit.utils.shell import run

log = get_logger("lnckit.deps")


@dataclass
class ToolSpec:
    """Declarative description of an external dependency."""
    name: str                       # canonical key, e.g. "hisat2"
    executables: tuple              # candidate binary names to search on PATH
    version_args: tuple = ("--version",)
    version_regex: str = r"(\d+\.\d+(?:\.\d+)?)"
    min_version: Optional[str] = None


@dataclass
class ResolvedTool:
    name: str
    path: str
    version: Optional[str]


# Registry of known tools. Adding a new backend = adding one entry here.
KNOWN_TOOLS = {
    # --- SRA download (the pipeline's entry point) --------------------------
    "prefetch":     ToolSpec("prefetch", ("prefetch",), ("--version",)),
    "fastq-dump":   ToolSpec("fastq-dump", ("fastq-dump",), ("--version",)),
    "fasterq-dump": ToolSpec("fasterq-dump", ("fasterq-dump",), ("--version",)),

    "fastqc":     ToolSpec("fastqc", ("fastqc",), ("--version",)),
    "multiqc":    ToolSpec("multiqc", ("multiqc",), ("--version",)),
    "fastp":      ToolSpec("fastp", ("fastp",), ("--version",)),
    "trimmomatic":ToolSpec("trimmomatic", ("trimmomatic",), ("--version",)),
    "cutadapt":   ToolSpec("cutadapt", ("cutadapt",), ("--version",)),
    "hisat2":     ToolSpec("hisat2", ("hisat2",), ("--version",), min_version="2.1.0"),
    "hisat2-build": ToolSpec("hisat2-build", ("hisat2-build",), ("--version",)),
    "star":       ToolSpec("star", ("STAR",), ("--version",)),
    "bowtie2":    ToolSpec("bowtie2", ("bowtie2",), ("--version",)),
    "bowtie2-build": ToolSpec("bowtie2-build", ("bowtie2-build",), ("--version",)),
    "bwa":        ToolSpec("bwa", ("bwa",), (), version_regex=r"Version:\s*(\d+\.\d+\.\d+)"),
    "samtools":   ToolSpec("samtools", ("samtools",), ("--version",), min_version="1.10"),
    "stringtie":  ToolSpec("stringtie", ("stringtie",), ("--version",), min_version="2.0.0"),
    "cufflinks":  ToolSpec("cufflinks", ("cufflinks",), ()),
    "gffread":    ToolSpec("gffread", ("gffread",), ("--version",)),
    "gffcompare": ToolSpec("gffcompare", ("gffcompare",), ("--version",)),
    "cpat":       ToolSpec("cpat", ("cpat.py", "cpat"), ("--version",)),
    "cpc2":       ToolSpec("cpc2", ("CPC2.py", "CPC2"), ()),
    "cnci":       ToolSpec("cnci", ("CNCI.py", "CNCI"), ()),
    "plek":       ToolSpec("plek", ("PLEK", "PLEK.py"), ()),
    "lncfinder":  ToolSpec("lncfinder", ("Rscript",), ("--version",),
                           version_regex=r"R scripting front-end version (\d+\.\d+\.\d+)"),
    "hmmscan":    ToolSpec("hmmscan", ("hmmscan",), ("-h",)),   # Pfam coding backend
    "hmmpress":   ToolSpec("hmmpress", ("hmmpress",), ("-h",)),
    "java":       ToolSpec("java", ("java",), ("-version",)),   # Weka (RNAplonc backend)
    "diamond":    ToolSpec("diamond", ("diamond",), ("version",), min_version="2.0.0"),
    "blastx":     ToolSpec("blastx", ("blastx",), ("-version",)),
    "blastn":     ToolSpec("blastn", ("blastn",), ("-version",)),
    "makeblastdb":ToolSpec("makeblastdb", ("makeblastdb",), ("-version",)),
    "cmscan":     ToolSpec("cmscan", ("cmscan",), ("-h",)),   # Infernal
    "cmpress":    ToolSpec("cmpress", ("cmpress",), ("-h",)),
    "feelnc":     ToolSpec("feelnc", ("FEELnc_filter.pl",), ()),
    "getorf":     ToolSpec("getorf", ("getorf",), ("-version",)),           # EMBOSS
    "transdecoder": ToolSpec("transdecoder", ("TransDecoder.LongOrfs",), ()),
    "orffinder":  ToolSpec("orffinder", ("ORFfinder",), ()),               # NCBI
}


def _extract_version(text: str, regex: str) -> Optional[str]:
    m = re.search(regex, text or "")
    return m.group(1) if m else None


def _ensure_managed_bin_on_path() -> None:
    """Prepend lnckit's managed bin dir (direct-installed tools) to PATH once."""
    try:
        from lnckit.core.installer import ensure_managed_bin_on_path
        ensure_managed_bin_on_path()
    except Exception:  # pragma: no cover - never let PATH setup break resolution
        pass


def resolve_tool(key: str) -> ResolvedTool:
    """Locate a tool on PATH and read its version. Raises DependencyError if absent."""
    _ensure_managed_bin_on_path()
    spec = KNOWN_TOOLS.get(key)
    if spec is None:
        raise DependencyError(f"Unknown tool '{key}' (not in dependency registry).")

    path = None
    for exe in spec.executables:
        found = shutil.which(exe)
        if found:
            path = found
            break
    if path is None:
        raise DependencyError(
            f"Required tool '{key}' not found on PATH "
            f"(looked for: {', '.join(spec.executables)}). "
            f"Install it via conda/bioconda or use the provided Docker/Singularity image."
        )

    ver = None
    if spec.version_args:
        try:
            res = run([path, *spec.version_args], tool=key, check=False)
            # 126/127 mean the shell found the file but could not execute it - a
            # missing shared library, a wrong architecture, a bad interpreter line.
            # Such a tool is worse than absent: it would pass a naive presence check
            # and then kill the pipeline hours later, so fail now, loudly.
            if res.returncode in (126, 127):
                raise BrokenToolError(
                    f"{key} is present at {path} but cannot be executed "
                    f"(exit {res.returncode}). Its installation is broken - reinstall "
                    f"it, e.g. 'conda install -c bioconda {key}'.\n"
                    f"  {(res.stderr or '').strip().splitlines()[0] if res.stderr else ''}",
                    tool=key, path=path,
                )
            # Only trust the output of a probe that actually succeeded. A failing
            # binary prints an error whose text often contains digits -
            # "libcrypto.so.1.0.0" would otherwise be scraped as version "1.0.0".
            if res.returncode == 0:
                ver = _extract_version(res.stdout + res.stderr, spec.version_regex)
            else:
                log.debug("%s: version probe exited %d; version unknown",
                          key, res.returncode)
        except BrokenToolError:
            raise
        except Exception:  # pragma: no cover - version probing is best-effort
            ver = None

    if spec.min_version and ver:
        try:
            if _v.parse(ver) < _v.parse(spec.min_version):
                raise OutdatedToolError(
                    f"{key} {ver} is installed at {path}, but lnckit requires "
                    f">= {spec.min_version}. Upgrade it (e.g. conda install -c bioconda "
                    f"'{key}>={spec.min_version}').",
                    tool=key, path=path, version=ver,
                )
        except _v.InvalidVersion:  # pragma: no cover
            pass

    log.debug("resolved %s -> %s (v%s)", key, path, ver)
    return ResolvedTool(key, path, ver)


def check_tools(keys) -> dict:
    """Resolve several tools, aggregating all missing ones into a single error."""
    resolved, missing = {}, []
    for k in keys:
        try:
            resolved[k] = resolve_tool(k)
        except DependencyError as e:
            missing.append(str(e))
    if missing:
        raise DependencyError("Missing/incompatible dependencies:\n  - " + "\n  - ".join(missing))
    return resolved
