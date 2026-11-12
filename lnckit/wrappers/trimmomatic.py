"""Trimmomatic + Cutadapt wrappers - alternative trimming backends to fastp.

Both are exposed behind the same `trim` stage interface (in1/out1[/in2/out2]) so
`tools.trim` can be switched in YAML without touching the pipeline.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register

# Trimmomatic's own recommended defaults; overridable via tool_params.trimmomatic.steps
DEFAULT_STEPS = ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]


@register("trim", "trimmomatic")
class TrimmomaticWrapper(BaseWrapper):
    """Trimmomatic PE/SE.

    `tool_params.trimmomatic` keys:
      adapters      : path to the adapter FASTA (e.g. TruSeq3-PE.fa). NOTE the
                      spelling - TruSeq, not TrueSeq; a wrong path makes
                      ILLUMINACLIP a no-op or a hard error.
      illuminaclip  : the ILLUMINACLIP suffix, default "2:30:10"
      steps         : list of trimming steps replacing DEFAULT_STEPS
    """

    requires = ("trimmomatic",)
    tool_name = "trimmomatic"

    def env(self, **_) -> Optional[dict]:
        """Give the JVM a real heap so Trimmomatic doesn't OOM.

        Trimmomatic's launcher honours `_JAVA_OPTIONS`. Without an explicit `-Xmx` the
        default heap is small and the parallel gzip compressor throws OutOfMemoryError on
        multi-threaded paired-end runs. Heap = `tool_params.trimmomatic.max_heap_gb` if
        set, else ~75% of the run's memory (leaving headroom), floor 2 GB.
        """
        import os
        heap = self.extra_params.get("max_heap_gb")
        if not heap:
            heap = max(2, int(int(self.extra_params.get("mem_gb", 8)) * 0.75))
        prev = os.environ.get("_JAVA_OPTIONS", "")
        return {"_JAVA_OPTIONS": (f"{prev} -Xmx{heap}g").strip()}

    def _unpaired(self, out: Path) -> Path:
        return out.parent / out.name.replace(".fastq.gz", ".unpaired.fastq.gz")

    @staticmethod
    def _bundled_adapter(paired: bool) -> Optional[str]:
        """Locate Trimmomatic's bundled TruSeq3-PE/SE.fa adapter (conda ships them under
        share/trimmomatic*/adapters/), so ILLUMINACLIP works out of the box."""
        import glob
        import shutil
        name = "TruSeq3-PE.fa" if paired else "TruSeq3-SE.fa"
        exe = shutil.which("trimmomatic")
        roots = []
        if exe:
            # <env>/bin/trimmomatic -> <env>/share/trimmomatic*/adapters
            roots.append(str(Path(exe).resolve().parents[1] / "share"))
        for env in ("CONDA_PREFIX",):
            import os
            if os.environ.get(env):
                roots.append(str(Path(os.environ[env]) / "share"))
        for root in roots:
            hits = glob.glob(f"{root}/trimmomatic*/adapters/{name}") + \
                glob.glob(f"{root}/trimmomatic*/adapters/{name.lower()}")
            if hits:
                return hits[0]
        return None

    def build_command(self, *, in1: Path, out1: Path,
                      in2: Optional[Path] = None, out2: Optional[Path] = None,
                      **_) -> List[str]:
        p = self.extra_params
        paired = bool(in2 and out2)
        cmd = ["trimmomatic", "PE" if paired else "SE", "-threads", self.threads]
        # Default -phred33 (modern Illumina). Override with tool_params.trimmomatic.phred:
        # 64 (old data) or 'none' to omit.
        phred = str(p.get("phred", "33"))
        if phred and phred != "none":
            cmd.append(f"-phred{phred}")

        if paired:
            cmd += [str(in1), str(in2), str(out1), str(self._unpaired(out1)),
                    str(out2), str(self._unpaired(out2))]
        else:
            cmd += [str(in1), str(out1)]

        # ILLUMINACLIP (adapter clipping): use the adapter the user set, else auto-find
        # Trimmomatic's bundled TruSeq3 file so adapters are ALWAYS clipped by default.
        # Disable entirely with tool_params.trimmomatic.adapters: none.
        adapters = p.get("adapters")
        if adapters == "none":
            adapters = None
        elif adapters:
            if not Path(adapters).exists():
                raise FileNotFoundError(
                    f"Trimmomatic adapter file not found: {adapters}. It ships with "
                    f"Trimmomatic under adapters/ (TruSeq3-PE.fa, TruSeq3-SE.fa, ...)."
                )
        else:
            adapters = self._bundled_adapter(paired)   # auto: TruSeq3-PE/SE.fa
        if adapters:
            # ILLUMINACLIP goes FIRST (Trimmomatic requires it before quality steps).
            # Standard suffix seedMismatches:palindromeClipThreshold:simpleClipThreshold.
            cmd.append(f"ILLUMINACLIP:{adapters}:{p.get('illuminaclip', '2:30:10')}")
        cmd += list(p.get("steps", DEFAULT_STEPS))
        return cmd

    def expected_outputs(self, *, out1: Path, out2: Optional[Path] = None, **_) -> List[Path]:
        return [Path(out1)] + ([Path(out2)] if out2 else [])


@register("trim", "cutadapt")
class CutadaptWrapper(BaseWrapper):
    """Cutadapt.

    `tool_params.cutadapt` keys: adapter_fwd (-a), adapter_rev (-A),
    quality_cutoff (-q, default 20), min_length (-m, default 36).
    """

    requires = ("cutadapt",)
    tool_name = "cutadapt"

    def build_command(self, *, in1: Path, out1: Path,
                      in2: Optional[Path] = None, out2: Optional[Path] = None,
                      **_) -> List[str]:
        p = self.extra_params
        cmd = ["cutadapt", "-j", self.threads,
               "-q", p.get("quality_cutoff", 20),
               "-m", p.get("min_length", 36)]
        if p.get("adapter_fwd"):
            cmd += ["-a", p["adapter_fwd"]]
        if in2 and p.get("adapter_rev"):
            cmd += ["-A", p["adapter_rev"]]
        cmd += ["-o", str(out1)]
        if in2 and out2:
            cmd += ["-p", str(out2), str(in1), str(in2)]
        else:
            cmd += [str(in1)]
        return cmd

    def expected_outputs(self, *, out1: Path, out2: Optional[Path] = None, **_) -> List[Path]:
        return [Path(out1)] + ([Path(out2)] if out2 else [])
