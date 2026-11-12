"""SRA Toolkit wrappers - turn public accessions into local FASTQ files.

Most published plant RNA-seq re-analyses start from an SRA/ENA accession list
rather than from local FASTQ, so accession -> FASTQ is treated as a first-class
pipeline stage rather than a manual pre-step.

Backends
--------
* ``prefetch``      - fetch the ``.sra`` archive into a local cache first. Strongly
  recommended: dumping straight from NCBI over the network is slow and fails
  mid-run on large runs, whereas a prefetched archive makes the dump local,
  restartable and reproducible.
* ``fastq-dump``    - the classic dumper. Slower, but writes ``.gz`` directly and
  is available in every SRA Toolkit build.
* ``fasterq-dump``  - multi-threaded, several times faster, but cannot compress;
  the fetch module gzips its output afterwards.

Split modes
-----------
``--split-files`` (default) writes ``_1``/``_2`` for properly paired spots and sends
orphaned mates to a separate unpaired file, so the paired files stay in strict
lock-step - which is what every downstream aligner assumes. ``--split-files``
(as used in many published command logs) is also supported, but it silently emits
mate files that can fall out of sync when a spot has a single read.
"""
from __future__ import annotations
import glob
import os
import re
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import List
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register

SPLIT_FLAGS = {
    "split-3": "--split-3",
    "split-files": "--split-files",
    "concatenate": "--concatenate-reads",
}


@lru_cache(maxsize=None)
def _sra_exe(name: str) -> str:
    """Resolve an SRA-Toolkit binary, PREFERRING a modern (>=3) build.

    Old 2.x builds (e.g. the conda `sra-tools=2.8.0`) reject `--output-directory`
    ('Unknown argument') and fail NCBI's TLS certificate check ('X509 ... not
    signed by trusted CA'), so accession download breaks. We therefore look past
    whatever is first on PATH: the package's own direct-download toolkit (in
    managed_bin_dir(), a known-good 3.x) is preferred, then any 3.x on PATH or in a
    conda env; a 2.x is used only if nothing newer exists.
    """
    def major(path: str) -> int:
        try:
            r = subprocess.run([path, "--version"], capture_output=True,
                               text=True, timeout=10)
            m = re.search(r"(\d+)\.(\d+)\.(\d+)", (r.stdout or "") + (r.stderr or ""))
            return int(m.group(1)) if m else 0
        except Exception:
            return -1

    cands: List[str] = []
    try:
        from lnckit.core.installer import managed_bin_dir
        mb = managed_bin_dir() / name              # package's own 3.x toolkit first
        if mb.exists():
            cands.append(str(mb))
    except Exception:
        pass
    w = shutil.which(name)
    if w and w not in cands:
        cands.append(w)
    home = os.path.expanduser("~")
    for pat in (f"{home}/.conda/envs/*/bin/{name}", f"/opt/conda/envs/*/bin/{name}"):
        for p in sorted(glob.glob(pat)):
            if p not in cands:
                cands.append(p)
    if not cands:
        return name
    best, best_mj = cands[0], major(cands[0])
    for c in cands:
        mj = major(c)
        if mj >= 3:                                # first modern build wins
            return c
        if mj > best_mj:
            best, best_mj = c, mj
    return best


@register("sra", "prefetch")
class PrefetchWrapper(BaseWrapper):
    """Download the .sra archive for one accession into a local cache."""

    requires = ("prefetch",)
    tool_name = "prefetch"

    def build_command(self, *, accession: str, cache_dir: Path,
                      max_size: str = "100G", **_) -> List[str]:
        # `-O` (not `--output-directory`) works on both old and new builds; the binary
        # is resolved to a modern (>=3) toolkit by _sra_exe (old 2.x can't download).
        return [_sra_exe("prefetch"), "-O", str(cache_dir),
                "--max-size", max_size, accession]

    def expected_outputs(self, *, accession: str, cache_dir: Path, **_) -> List[Path]:
        # prefetch writes <cache>/<ACC>/<ACC>.sra (older builds: <cache>/<ACC>.sra);
        # the fetch module resolves whichever exists, so nothing is asserted here.
        return []


@register("sra", "fastq-dump")
class FastqDumpWrapper(BaseWrapper):
    """Classic dumper: writes gzipped FASTQ directly."""

    requires = ("fastq-dump",)
    tool_name = "fastq-dump"

    def build_command(self, *, accession: str, outdir: Path,
                      split: str = "split-files", gzip: bool = True,
                      sra_file: Path = None, **_) -> List[str]:
        cmd = [_sra_exe("fastq-dump"), SPLIT_FLAGS.get(split, "--split-files"),
               "-O", str(outdir)]
        if gzip:
            cmd.append("--gzip")
        cmd += [str(sra_file) if sra_file else accession]
        return cmd

    def expected_outputs(self, **_) -> List[Path]:
        # The number of FASTQ files depends on the library layout (single vs
        # paired), which is only known after the dump; modules/fetch.py globs and
        # validates the result instead.
        return []


@register("sra", "fasterq-dump")
class FasterqDumpWrapper(BaseWrapper):
    """Multi-threaded dumper (uncompressed output; gzipped by the fetch module)."""

    requires = ("fasterq-dump",)
    tool_name = "fasterq-dump"

    def build_command(self, *, accession: str, outdir: Path, tmp_dir: Path = None,
                      split: str = "split-files", sra_file: Path = None, **_) -> List[str]:
        cmd = [_sra_exe("fasterq-dump"), SPLIT_FLAGS.get(split, "--split-files"),
               "--threads", self.threads, "-O", str(outdir)]
        if tmp_dir:
            cmd += ["--temp", str(tmp_dir)]
        cmd += [str(sra_file) if sra_file else accession]
        return cmd

    def expected_outputs(self, **_) -> List[Path]:
        return []
