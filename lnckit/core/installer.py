"""Dependency doctor + auto-installer for lnckit's external tools.

lnckit wraps established bioinformatics software instead of reimplementing it,
so a fresh machine will usually be missing several executables. This module can

  1. *report* exactly which required tools are present or missing (``doctor``), and
  2. *install* the missing ones, choosing the best available method per tool.

Install methods, in order of preference
---------------------------------------
* **conda / mamba (bioconda)** - the correct, reproducible way to obtain the whole
  stack (HISAT2, STAR, samtools, EMBOSS/getorf, FEELnc, Infernal, ...). Used first
  whenever a conda/mamba executable is on PATH.
* **pip** - for the handful of tools published as Python packages (CPAT, cutadapt).
* **direct** - download an official, self-contained release binary from the tool's
  GitHub releases and drop it into a lnckit-managed ``bin`` directory. Only used
  for tools that genuinely ship a single static binary (DIAMOND, StringTie, STAR,
  gffread, gffcompare, fastp). This is the "download it directly" fallback for
  machines without conda; it is deliberately *not* attempted for tools that need
  compilation or a runtime (Perl/Java/EMBOSS suite), where conda is the only sane
  route.

Anything installed by the ``direct`` method lands in ``managed_bin_dir()``, which
lnckit automatically prepends to ``PATH`` at runtime, so subsequent runs find it.
"""
from __future__ import annotations

import json
import os
import platform
import shutil
import stat
import sys
import tarfile
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from lnckit.core.dependencies import KNOWN_TOOLS, resolve_tool
from lnckit.core.exceptions import (
    BrokenToolError, DependencyError, OutdatedToolError, LncKitError)
from lnckit.core.logging_utils import get_logger

log = get_logger("lnckit.installer")


# --------------------------------------------------------------------------- #
# Managed location for directly-downloaded binaries
# --------------------------------------------------------------------------- #
def managed_bin_dir() -> Path:
    """Directory where ``direct`` installs place binaries (XDG-friendly)."""
    base = os.environ.get("LNCKIT_HOME") or os.environ.get("XDG_DATA_HOME")
    root = Path(base) if base else Path.home() / ".local" / "share"
    d = root / "lnckit" / "bin"
    return d


def ensure_managed_bin_on_path() -> None:
    """Prepend the managed bin dir to PATH for this process (idempotent)."""
    d = str(managed_bin_dir())
    parts = os.environ.get("PATH", "").split(os.pathsep)
    if os.path.isdir(d) and d not in parts:
        os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")


def managed_src_dir() -> Path:
    """Directory where ``source`` installs clone/unpack their trees."""
    base = os.environ.get("LNCKIT_HOME") or os.environ.get("XDG_DATA_HOME")
    root = Path(base) if base else Path.home() / ".local" / "share"
    return root / "lnckit" / "src"


# --------------------------------------------------------------------------- #
# Platform detection
# --------------------------------------------------------------------------- #
def _os_key() -> str:
    s = platform.system().lower()
    if s == "linux":
        return "linux"
    if s == "darwin":
        return "mac"
    return s


def _is_linux() -> bool:
    return _os_key() == "linux"


# --------------------------------------------------------------------------- #
# Install recipes
# --------------------------------------------------------------------------- #
@dataclass
class DirectRecipe:
    """How to fetch a self-contained binary from a GitHub release."""
    repo: str                       # "owner/name"
    asset_regex_linux: str          # regex matched against release asset names
    asset_regex_mac: Optional[str] = None
    binaries: tuple = ()            # basenames to extract from the archive
    # If the asset is itself the raw binary (not an archive), set raw_binary_name:
    raw_binary_name: Optional[str] = None


@dataclass
class SourceRecipe:
    """How to install a script/source tool that isn't a single prebuilt binary.

    Used for tools distributed as a plain tarball (cufflinks), a git repo with a
    build step (CPC2), or a SourceForge archive that compiles on install (PLEK).
    The tool is fetched into managed_src_dir()/<tool>, optionally built, and the
    executables in ``exposes`` are symlinked into managed_bin_dir() (which is on
    PATH), so they resolve exactly like conda/direct installs.
    """
    url: Optional[str] = None       # tarball/zip URL (mutually exclusive with git)
    git: Optional[str] = None       # git repo to clone (shallow)
    build: Optional[str] = None     # shell command run in the fetched source dir
    exposes: tuple = ()             # executable basenames to find + put on PATH
    aliases: dict = field(default_factory=dict)  # extra_name -> exposed_name symlink
    interpreter: Optional[str] = None  # conda env name whose python runs the script
                                       # (for Python-2 tools like CPC2/PLEK)
    interpreter_spec: tuple = ()       # conda packages to build that env if missing


@dataclass
class InstallRecipe:
    tool: str                       # canonical KNOWN_TOOLS key
    conda: Optional[str] = None     # bioconda/conda-forge package spec
    pip: Optional[str] = None       # PyPI package name
    direct: Optional[DirectRecipe] = None
    source: Optional[SourceRecipe] = None
    system: Optional[str] = None    # OS package name for sudo <pkg-mgr> install
    note: str = ""

    def methods(self) -> list:
        m = []
        if self.conda:
            m.append("conda")
        if self.pip:
            m.append("pip")
        if self.source:
            m.append("source")
        if self.direct:
            m.append("direct")
        if self.system:
            m.append("system")
        return m


# SRA Toolkit: one download provides all three executables. We use the official
# NCBI tarball (self-contained static binaries) rather than conda on purpose: the
# bioconda sra-tools build pins an ancient htslib that downgrades samtools in the
# shared env. Builds older than 3.0 also fail NCBI's current TLS certs
# ("X509 - Certificate verification failed"), which the 3.x tarball avoids.
_SRA_URL = ("https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/"
            "sratoolkit.3.4.1-ubuntu64.tar.gz")
_SRA_BINS = ("prefetch", "fastq-dump", "fasterq-dump")

# Curated recipe table. conda package names follow bioconda.
INSTALL_RECIPES: dict[str, InstallRecipe] = {
    "prefetch":    InstallRecipe(
        "prefetch",
        source=SourceRecipe(url=_SRA_URL, exposes=_SRA_BINS),
        note="official NCBI SRA Toolkit tarball; also provides fastq-dump"),
    "fastq-dump":  InstallRecipe(
        "fastq-dump",
        source=SourceRecipe(url=_SRA_URL, exposes=_SRA_BINS)),
    "fasterq-dump":InstallRecipe(
        "fasterq-dump",
        source=SourceRecipe(url=_SRA_URL, exposes=_SRA_BINS)),

    "fastqc":      InstallRecipe("fastqc", conda="fastqc"),
    "multiqc":     InstallRecipe("multiqc", conda="multiqc", pip="multiqc"),
    "fastp":       InstallRecipe(
        "fastp", conda="fastp", pip=None,
        direct=DirectRecipe("OpenGene/fastp",
                            asset_regex_linux=r"^fastp$|^fastp\.gz$|^fastp\.\d",
                            binaries=("fastp",), raw_binary_name="fastp"),
        note="direct download needs a prebuilt asset; otherwise use conda",
    ),
    "trimmomatic": InstallRecipe("trimmomatic", conda="trimmomatic"),
    "cutadapt":    InstallRecipe("cutadapt", conda="cutadapt", pip="cutadapt"),
    "hisat2":      InstallRecipe("hisat2", conda="hisat2",
                                 note="also provides hisat2-build"),
    "hisat2-build":InstallRecipe("hisat2-build", conda="hisat2"),
    "star":        InstallRecipe(
        "star", conda="star",
        direct=DirectRecipe("alexdobin/STAR",
                            asset_regex_linux=r"STAR.*Linux.*x86_64.*\.zip|Linux_x86_64",
                            binaries=("STAR",)),
    ),
    "cufflinks":   InstallRecipe(
        "cufflinks",
        source=SourceRecipe(
            url="http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/"
                "cufflinks-2.2.1.Linux_x86_64.tar.gz",
            exposes=("cufflinks", "cuffcompare", "gtf_to_sam"),
        ),
        note="prebuilt Linux x86_64 tarball; no build step",
    ),
    "bowtie2":     InstallRecipe("bowtie2", conda="bowtie2"),
    "bowtie2-build":InstallRecipe("bowtie2-build", conda="bowtie2"),
    "bwa":         InstallRecipe("bwa", conda="bwa"),
    "samtools":    InstallRecipe(
        "samtools", conda="samtools>=1.15",
        source=SourceRecipe(
            url="https://github.com/samtools/samtools/releases/download/"
                "1.21/samtools-1.21.tar.bz2",
            # release tarball bundles htslib; needs dev headers (see SYSTEM_BUILD_DEPS)
            build="bash -c './configure && make'",
            exposes=("samtools",),
        ),
        system="samtools",
        note="conda preferred; source build auto-pulls build deps via the system "
             "package manager, or install directly with 'sudo dnf install samtools'"),
    "stringtie":   InstallRecipe(
        "stringtie", conda="stringtie",
        direct=DirectRecipe("gpertea/stringtie",
                            asset_regex_linux=r".*Linux_x86_64.*\.tar\.gz$",
                            binaries=("stringtie",)),
    ),
    "gffread":     InstallRecipe(
        "gffread", conda="gffread",
        direct=DirectRecipe("gpertea/gffread",
                            asset_regex_linux=r".*Linux_x86_64.*\.tar\.gz$",
                            binaries=("gffread",)),
    ),
    "gffcompare":  InstallRecipe(
        "gffcompare", conda="gffcompare",
        direct=DirectRecipe("gpertea/gffcompare",
                            asset_regex_linux=r".*Linux_x86_64.*\.tar\.gz$",
                            binaries=("gffcompare",)),
    ),
    "hmmsearch":   InstallRecipe("hmmsearch", conda="hmmer"), # Pfam domain scan (fast path)
    "hmmscan":     InstallRecipe("hmmscan", conda="hmmer"),   # Pfam domain scan
    "hmmpress":    InstallRecipe("hmmpress", conda="hmmer"),
    "java":        InstallRecipe("java", conda="openjdk"),    # Weka (RNAplonc backend)
    # pip only (NOT conda): bioconda's solver can pull the old CPAT 1.2.x, which is
    # Python-2 source. Installed into the Python-3 lnckit env its shebang runs it under
    # python3 and it crashes on py2 'print >>sys.stderr' syntax. The PyPI 'CPAT' package
    # is the maintained Python-3 v3 (same -x/-d/-g/-o CLI; parse_cpat reads its output).
    "cpat":        InstallRecipe("cpat", pip="CPAT"),
    "cpc2":        InstallRecipe(
        "cpc2",
        source=SourceRecipe(
            git="https://github.com/gao-lab/CPC2_standalone.git",
            # CPC2 ships libsvm as a tarball; unpack it, then compile once.
            build="bash -c 'cd libs/libsvm && tar xzf libsvm-*.tar.gz "
                  "&& cd libsvm-*/ && make clean && make'",
            exposes=("CPC2.py",),
            interpreter="lnckit_py2",   # CPC2 is a Python-2 tool
            interpreter_spec=("python=2.7", "biopython=1.76", "numpy=1.16"),
        ),
        note="git clone + libsvm build; runs under the lnckit_py2 env",
    ),
    "plek":        InstallRecipe(
        "plek",
        source=SourceRecipe(
            url="https://downloads.sourceforge.net/project/plek/PLEK.1.2.tar.gz",
            # PLEK compiles its SVM helpers on install (Python 2 tool).
            build="bash -c 'cd PLEK.1.2 && python PLEK_setup.py || true'",
            exposes=("PLEK.py",),
            aliases={"PLEK": "PLEK.py"},
            interpreter="lnckit_py2",   # PLEK is a Python-2 tool
            interpreter_spec=("python=2.7", "biopython=1.76", "numpy=1.16"),
        ),
        note="SourceForge tarball; runs under the lnckit_py2 env",
    ),
    "cnci":        InstallRecipe(
        "cnci",
        source=SourceRecipe(
            git="https://github.com/www-bioinfo-org/CNCI.git",
            # CNCI ships libsvm-3.0 as source in a zip; unpack AND compile it so the
            # svm-scale / svm-predict binaries CNCI shells out to exist.
            build="bash -c 'unzip -o libsvm-3.0.zip && make -C libsvm-3.0'",
            exposes=("CNCI.py",),
            interpreter="lnckit_py2",   # CNCI is a Python-2 tool
            interpreter_spec=("python=2.7", "biopython=1.76", "numpy=1.16"),
        ),
        note="git clone + unzip libsvm; runs under the lnckit_py2 env",
    ),
    "lncfinder":   InstallRecipe(
        "lncfinder", conda="r-base",
        note="LncFinder is CRAN-only (no conda package). setup.sh builds the isolated "
             "lnckit_lncfinder env: conda for the heavy deps (r-caret/r-seqinr/r-e1071) "
             "then install.packages('LncFinder') from CRAN. The wrapper auto-finds it."),
    "diamond":     InstallRecipe(
        "diamond", conda="diamond",
        direct=DirectRecipe("bbuchfink/diamond",
                            asset_regex_linux=r"diamond-linux64\.tar\.gz$",
                            binaries=("diamond",)),
    ),
    "blastx":      InstallRecipe("blastx", conda="blast"),
    "blastn":      InstallRecipe("blastn", conda="blast"),
    "makeblastdb": InstallRecipe("makeblastdb", conda="blast"),
    "cmscan":      InstallRecipe("cmscan", conda="infernal"),
    "cmpress":     InstallRecipe("cmpress", conda="infernal"),
    "feelnc":      InstallRecipe("feelnc", conda="feelnc",
                                 note="Perl; conda only"),
    "getorf":      InstallRecipe("getorf", conda="emboss",
                                 note="EMBOSS suite; conda only"),
    "transdecoder":InstallRecipe("transdecoder", conda="transdecoder"),
    "orffinder":   InstallRecipe("orffinder", conda="orffinder",
                                 note="NCBI ORFfinder"),
}


# --------------------------------------------------------------------------- #
# Which tools does a given configuration actually require?
# --------------------------------------------------------------------------- #
def required_tools(cfg) -> list:
    """Map the selected backends in a Config to canonical tool keys."""
    keys: list = []

    def add(k):
        if k and k not in keys:
            keys.append(k)

    tools = getattr(cfg, "tools", {}) or {}

    # SRA download stage - only required when the run starts from accessions.
    try:
        accs = cfg.get("inputs.sra_accessions") or cfg.get("inputs.sra_accession_file")
    except Exception:  # pragma: no cover - defensive
        accs = None
    if accs:
        sra = cfg.get("sra", {}) or {}
        add(sra.get("tool", "fastq-dump"))
        if sra.get("prefetch", True):
            add("prefetch")

    aligner = tools.get("aligner")
    add(aligner)
    if aligner == "hisat2":
        add("hisat2-build")
    elif aligner == "bowtie2":
        add("bowtie2-build")
    add("samtools")
    add(tools.get("qc"))
    add(tools.get("trim"))
    add(tools.get("assembler"))
    add("stringtie")          # stringtie --merge is the merge backend for both assemblers
    add("gffread")
    add("gffcompare")
    # coding_potential may be a single tool or an ensemble list (CPAT + LncFinder ...)
    cp = tools.get("coding_potential")
    for c in (cp if isinstance(cp, (list, tuple)) else [cp]):
        if c == "rnaplonc":
            add("java")           # vendored REPTree/Weka classifier needs a JVM
        else:
            add(c)
    add(tools.get("protein_filter"))
    if tools.get("protein_filter") == "blastx":
        add("makeblastdb")    # the protein DB must be formatted before blastx runs
    housekeeping = tools.get("housekeeping")
    add(housekeeping)
    if housekeeping == "infernal":
        add("cmscan")
        add("cmpress")        # Rfam.cm must be pressed before cmscan can use it
    elif housekeeping == "blastn":
        add("makeblastdb")
    add(tools.get("classifier"))

    # ORF backend from filters.orf.tool (only if the ORF filter is enabled)
    try:
        orf = cfg.get("filters.orf", {}) or {}
        if orf.get("enabled", True):
            add(orf.get("tool", "getorf"))
    except Exception:  # pragma: no cover - defensive
        pass

    # keep only tools we know how to check
    return [k for k in keys if k in KNOWN_TOOLS]


# --------------------------------------------------------------------------- #
# Doctor: present/missing report
# --------------------------------------------------------------------------- #
@dataclass
class ToolStatus:
    tool: str
    present: bool
    version: Optional[str]
    path: Optional[str]
    methods: list = field(default_factory=list)  # available install methods
    #: "present" | "missing" | "outdated" - an outdated tool is on PATH but unusable,
    #: and needs upgrading rather than installing.
    state: str = "missing"
    detail: str = ""


def doctor(cfg=None, tools: Optional[list] = None) -> list:
    """Return status for each required (or given) tool. Makes managed bin visible."""
    ensure_managed_bin_on_path()
    if tools is None:
        tools = required_tools(cfg) if cfg is not None else list(KNOWN_TOOLS)
    out = []
    for key in tools:
        recipe = INSTALL_RECIPES.get(key)
        methods = recipe.methods() if recipe else []
        try:
            rt = resolve_tool(key)
            out.append(ToolStatus(key, True, rt.version, rt.path, methods,
                                  state="present"))
        except OutdatedToolError as e:
            out.append(ToolStatus(key, False, e.version, e.path, methods,
                                  state="outdated", detail=str(e)))
        except BrokenToolError as e:
            out.append(ToolStatus(key, False, None, e.path, methods,
                                  state="broken", detail=str(e)))
        except DependencyError as e:
            out.append(ToolStatus(key, False, None, None, methods,
                                  state="missing", detail=str(e)))
    return out


# --------------------------------------------------------------------------- #
# Installer back-ends
# --------------------------------------------------------------------------- #
def conda_exe() -> Optional[str]:
    """Prefer mamba (faster solver) then conda/micromamba."""
    for exe in ("mamba", "micromamba", "conda"):
        p = shutil.which(exe)
        if p:
            return p
    return None


def conda_env_python(env_name: str) -> Optional[Path]:
    """Locate the `python` of a named conda env, or None if it can't be found.

    Named envs do NOT always live under <base>/envs: on a shared/read-only base
    (typical HPC, e.g. /cm/shared/apps/anaconda3) conda puts a user's named envs in
    ~/.conda/envs instead, and a user may set CONDA_ENVS_DIRS. We check all of these.
    """
    exe = conda_exe()
    env_dirs: list[Path] = []

    # 0. the envs dir that holds the RUNNING lnckit env - its companion env (lnckit_py2)
    #    is built right beside it, so this is correct for conda AND micromamba installs
    #    (sys.prefix is .../envs/lnckit, its parent is the envs dir).
    env_dirs.append(Path(sys.prefix).parent)

    # 1. explicit envs directories (these ARE the parent of <name>/bin/python)
    for var in ("CONDA_ENVS_PATH", "CONDA_ENVS_DIRS"):
        if os.environ.get(var):
            env_dirs += [Path(p) for p in os.environ[var].split(os.pathsep) if p]
    # 2. ~/.conda/envs - where a read-only base lands a user's named envs
    env_dirs.append(Path.home() / ".conda" / "envs")

    # 3. <base>/envs for every plausible conda/mamba base
    bases: list[Path] = []
    if exe:
        # <base>/bin/conda or <base>/condabin/conda -> base is two levels up
        bases.append(Path(exe).resolve().parent.parent)
    for env in ("CONDA_PREFIX", "CONDA_PREFIX_1", "MAMBA_ROOT_PREFIX"):
        if os.environ.get(env):
            bases.append(Path(os.environ[env]))
    bases += [Path.home() / d for d in
              ("anaconda3", "miniconda3", "miniforge3", "micromamba")]
    bases.append(Path("/opt/conda"))
    env_dirs += [b / "envs" for b in bases]

    seen = set()
    for d in env_dirs:
        if d in seen:
            continue
        seen.add(d)
        cand = d / env_name / "bin" / "python"
        if cand.exists():
            return cand
    return None


def _env_manager() -> Optional[str]:
    """The manager to build a companion env with - prefer the fast one that owns the
    running env (micromamba/mamba) so the new env lands in the SAME envs dir; else conda."""
    for m in ("micromamba", "mamba", "conda"):
        p = shutil.which(m)
        if p:
            return p
    return None


def _create_companion_env(name: str, spec: tuple) -> Optional[Path]:
    """Create a companion env (e.g. lnckit_py2 for CPC2/CNCI/PLEK) right BESIDE the running
    lnckit env, with whichever manager is available, so its python is at the same envs dir
    as the rest - correct whether LncKit was installed with conda OR micromamba. A `-p`
    prefix env (not `-n`) pins the location exactly. Returns its python path, or None."""
    import subprocess
    mgr = _env_manager()
    if not (spec and mgr):
        return None
    prefix = Path(sys.prefix).parent / name          # sibling of the lnckit env
    log.info("creating %s env at %s (via %s): %s",
             name, prefix, os.path.basename(mgr), " ".join(spec))
    try:
        subprocess.run([mgr, "create", "-y", "-p", str(prefix),
                        "--override-channels", "-c", "conda-forge", "-c", "bioconda",
                        *spec], check=True)          # noqa: S603
    except Exception as e:                            # noqa: BLE001 - offline/solve failure
        log.warning("could not create %s env (%s)", name, e)
        return None
    py = prefix / "bin" / "python"
    return py if py.exists() else None


def _gh_headers() -> dict:
    h = {"User-Agent": "lnckit-installer", "Accept": "application/vnd.github+json"}
    tok = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if tok:
        h["Authorization"] = f"Bearer {tok}"
    return h


def _http_json(url: str) -> dict:
    import urllib.error
    req = urllib.request.Request(url, headers=_gh_headers())
    try:
        with urllib.request.urlopen(req, timeout=30) as r:  # noqa: S310 (trusted host)
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 403 and "rate limit" in (e.headers.get("x-ratelimit-remaining", "") + str(e.reason)).lower() \
                or e.code == 403:
            raise LncKitError(
                "GitHub API rate limit hit (unauthenticated = 60 requests/hour). "
                "Set a GITHUB_TOKEN environment variable to raise the limit, or "
                "install via conda instead (--method conda)."
            ) from e
        raise LncKitError(f"GitHub API error {e.code} for {url}: {e.reason}") from e


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "lnckit-installer"})
    log.info("downloading %s", url)
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:  # noqa: S310
        shutil.copyfileobj(r, f)
    return dest


def _github_asset_url(recipe: DirectRecipe) -> tuple:
    """Return (asset_name, download_url) for the latest matching release asset."""
    import re
    data = _http_json(f"https://api.github.com/repos/{recipe.repo}/releases/latest")
    pattern = recipe.asset_regex_linux if _is_linux() else (
        recipe.asset_regex_mac or recipe.asset_regex_linux)
    rx = re.compile(pattern)
    for asset in data.get("assets", []):
        if rx.search(asset["name"]):
            return asset["name"], asset["browser_download_url"]
    raise LncKitError(
        f"No release asset matching /{pattern}/ found for {recipe.repo} "
        f"(tag {data.get('tag_name')}). Use conda instead: "
        f"the tool may not ship a prebuilt {_os_key()} binary."
    )


def _extract_binaries(archive: Path, wanted: tuple, dest_dir: Path) -> list:
    """Extract the requested binary basenames from a tar/zip archive."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    installed = []
    wanted_set = set(wanted)

    def _emit(name: str, fileobj) -> None:
        base = os.path.basename(name)
        if base in wanted_set:
            target = dest_dir / base
            with open(target, "wb") as out:
                shutil.copyfileobj(fileobj, out)
            target.chmod(target.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
            installed.append(base)

    if tarfile.is_tarfile(archive):
        with tarfile.open(archive) as tf:
            for m in tf.getmembers():
                if m.isfile() and os.path.basename(m.name) in wanted_set:
                    f = tf.extractfile(m)
                    if f:
                        _emit(m.name, f)
    elif zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as zf:
            for n in zf.namelist():
                if os.path.basename(n) in wanted_set and not n.endswith("/"):
                    with zf.open(n) as f:
                        _emit(n, f)
    else:
        raise LncKitError(f"Unrecognised archive format: {archive.name}")

    missing = wanted_set - set(installed)
    if missing:
        raise LncKitError(
            f"Archive {archive.name} did not contain expected binaries: "
            f"{', '.join(sorted(missing))}")
    return installed


def install_direct(tool: str, dry_run: bool = False) -> str:
    """Download a self-contained release binary into the managed bin dir."""
    recipe = INSTALL_RECIPES.get(tool)
    if not recipe or not recipe.direct:
        raise LncKitError(f"No direct-download recipe for '{tool}'.")
    d = recipe.direct
    name, url = _github_asset_url(d)
    bindir = managed_bin_dir()
    if dry_run:
        return f"direct: would download {name} from {d.repo} -> {bindir}"

    bindir.mkdir(parents=True, exist_ok=True)
    tmp = bindir / name
    _download(url, tmp)
    try:
        # raw single-binary asset (e.g. fastp), possibly gzipped, vs. archive
        is_tar = tarfile.is_tarfile(tmp)
        if d.raw_binary_name and not (is_tar or zipfile.is_zipfile(tmp)):
            target = bindir / d.raw_binary_name
            if name.endswith(".gz"):
                import gzip
                with gzip.open(tmp, "rb") as fin, open(target, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                tmp.unlink(missing_ok=True)
            else:
                shutil.move(str(tmp), target)
            target.chmod(target.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
            installed = [d.raw_binary_name]
        else:
            installed = _extract_binaries(tmp, d.binaries, bindir)
            tmp.unlink(missing_ok=True)
    finally:
        if tmp.exists() and tmp.is_file() and tmp.suffix in {".gz", ".zip", ".tar", ".bz2", ".xz"}:
            tmp.unlink(missing_ok=True)
    ensure_managed_bin_on_path()
    return f"direct: installed {', '.join(installed)} -> {bindir}"


def _extract_all(archive: Path, dest_dir: Path) -> None:
    """Unpack a whole tar/zip archive into dest_dir (preserving its layout)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    if tarfile.is_tarfile(archive):
        with tarfile.open(archive) as tf:
            tf.extractall(dest_dir)  # noqa: S202 - trusted, curated recipe URLs
    elif zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(dest_dir)
    else:
        raise LncKitError(f"Unrecognised archive format: {archive.name}")


def _find_in_tree(root: Path, basename: str) -> Optional[Path]:
    """Return the shallowest file named `basename` under root, or None."""
    matches = sorted(root.rglob(basename), key=lambda p: len(p.parts))
    for m in matches:
        if m.is_file():
            return m
    return None


def _py2_shim_text(env_name: str, script: Path, resolved: Optional[Path]) -> str:
    """Build a launcher that resolves the interpreter env's python AT RUN TIME.

    A named env's location differs between conda and micromamba installs, and an
    ABSOLUTE interpreter path baked in at install time goes stale the moment the
    env is rebuilt, moved, or (on a shared HPC) lands in a different envs dir than
    when the shim was written - the tool then dies instantly with
    `bad interpreter: No such file or directory`. So instead of hardcoding one
    path, the shim probes the usual locations on every call and execs the first
    python that actually exists. An override env var is honoured first.
    """
    import re
    beside = Path(sys.prefix).parent / env_name / "bin" / "python"   # beside the lnckit env
    var = re.sub(r"[^A-Za-z0-9]", "_", env_name).upper() + "_PYTHON"   # e.g. LNCKIT_PY2_PYTHON
    cands = ['"${%s:-}"' % var]                                      # 1. explicit override
    if resolved:
        cands.append('"%s"' % resolved)                             # 2. what we resolved now
    cands += [
        '"%s"' % beside,                                            # 3. sibling of lnckit env
        '"${CONDA_PREFIX%%/*}/%s/bin/python"' % env_name,          # 4. .../envs/lnckit -> .../envs
        '"${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/%s/bin/python"' % env_name,
        '"$HOME/.conda/envs/%s/bin/python"' % env_name,
        '"$HOME/micromamba/envs/%s/bin/python"' % env_name,
    ]
    loop = " \\\n        ".join(cands)
    return (
        '#!/usr/bin/env bash\n'
        '# lnckit launcher: run this script under the "%s" env python, resolved AT RUN\n'
        '# TIME. Do NOT hardcode the path - it differs across conda/micromamba and a\n'
        '# baked path breaks if the env is rebuilt or moved (bad interpreter at run).\n'
        '_lk_script="%s"\n'
        'for _lk_py in \\\n        %s ; do\n'
        '    if [ -n "$_lk_py" ] && [ -x "$_lk_py" ]; then exec "$_lk_py" "$_lk_script" "$@"; fi\n'
        'done\n'
        'echo "lnckit: could not find the %s env python for $(basename "$_lk_script")." >&2\n'
        'echo "        The %s env is missing - re-run ./lnckit.sh to build it, or create it:" >&2\n'
        'echo "          <mgr> create -y -p <envs_dir>/%s --override-channels -c conda-forge -c bioconda python=2.7 biopython=1.76 numpy=1.16" >&2\n'
        'echo "        (or point %s=/path/to/%s/bin/python at an existing one)." >&2\n'
        'exit 127\n'
    ) % (env_name, script, loop, env_name, env_name, env_name, var, env_name)


def repair_py2_shims() -> list:
    """Heal any managed-bin Python-2 launcher whose baked interpreter path no longer
    exists, rewriting it as the run-time-resolving shim (_py2_shim_text).

    Older versions (and setup.sh's CPC2 block) wrote a fixed `exec "<python>" "<script>"
    "$@"` line. If that env later moves, is rebuilt, or lands in a different envs dir than
    when the shim was written (classic on a shared HPC: `~/.conda/envs/...` vs micromamba),
    the tool dies instantly with `bad interpreter: No such file`. The wrapped script is
    still on disk, so we just rewrite the launcher in place - NO rebuild, NO download.
    Returns the basenames that were repaired."""
    import re as _re
    bindir = managed_bin_dir()
    repaired: list = []
    if not bindir.is_dir():
        return repaired
    for p in sorted(bindir.iterdir()):
        try:
            if p.is_symlink() or not p.is_file():
                continue                       # aliases are symlinks; plain tools too
            with p.open("rb") as fh:
                head_bytes = fh.read(4096)     # shims are tiny; don't slurp real binaries
            if not head_bytes.startswith(b"#!/usr/bin/env bash"):
                continue
            head = head_bytes.decode("utf-8", "replace")
            if "_lk_script=" in head or "_lk_py" in head:
                continue                       # already a run-time-resolving shim
            m = _re.search(r'^exec\s+"([^"]+)"\s+"([^"]+)"\s+"\$@"\s*$', head, _re.M)
            if not m:
                continue                       # not an old single-exec py2 shim
            interp, script = m.group(1), m.group(2)
            if "$" in interp or Path(interp).exists():
                continue                       # interpreter is a var, or it still works
            # derive the env name from the dead path (.../envs/<name>/bin/python)
            parts = Path(interp).parts
            env_name = parts[-3] if len(parts) >= 3 and parts[-1] == "python" else "lnckit_py2"
            p.write_text(_py2_shim_text(env_name, Path(script), None))
            p.chmod(0o755)
            repaired.append(p.name)
        except Exception as e:                 # noqa: BLE001 - one bad file must not stop the rest
            log.warning("could not inspect/repair shim %s (%s)", getattr(p, "name", p), e)
    return repaired


def install_from_source(tool: str, dry_run: bool = False) -> str:
    """Fetch (tarball or git), optionally build, and expose a script tool.

    The tree lands in managed_src_dir()/<tool>; the executables named in the
    recipe's ``exposes`` are symlinked into managed_bin_dir() so they are found
    on PATH like any other install.
    """
    import subprocess
    recipe = INSTALL_RECIPES.get(tool)
    if not recipe or not recipe.source:
        raise LncKitError(f"No source recipe for '{tool}'.")
    s = recipe.source
    srcroot = managed_src_dir() / tool
    bindir = managed_bin_dir()
    if dry_run:
        how = f"git clone {s.git}" if s.git else f"download {s.url}"
        build = f", build ({s.build.split(chr(39))[0].strip()}...)" if s.build else ""
        return (f"source: would {how} -> {srcroot}{build}, "
                f"expose {', '.join(s.exposes)} -> {bindir}")

    srcroot.parent.mkdir(parents=True, exist_ok=True)
    if srcroot.exists():
        shutil.rmtree(srcroot)

    # 1. fetch
    if s.git:
        log.info("cloning %s", s.git)
        subprocess.run(["git", "clone", "--depth", "1", s.git, str(srcroot)],  # noqa: S603,S607
                       check=True)
    elif s.url:
        srcroot.mkdir(parents=True, exist_ok=True)
        tmp = srcroot / "_download"
        _download(s.url, tmp)
        _extract_all(tmp, srcroot)
        tmp.unlink(missing_ok=True)
    else:
        raise LncKitError(f"source recipe for '{tool}' has neither git nor url")

    # 2. build (optional). A compile may need system dev headers first.
    if s.build:
        if tool in SYSTEM_BUILD_DEPS:
            try:
                dep_line = install_build_deps(tool)
                if dep_line:
                    log.info("%s build deps: %s", tool, dep_line)
            except Exception as e:      # noqa: BLE001 - deps may already be present
                log.warning("could not install build deps for %s (%s); "
                            "continuing", tool, e)
        log.info("building %s: %s", tool, s.build)
        subprocess.run(s.build, shell=True, cwd=str(srcroot), check=True)  # noqa: S602

    # 3. expose executables on PATH. A plain tool is symlinked; a Python-2 tool
    #    (recipe.interpreter set) gets a shim that runs it under that env's python.
    bindir.mkdir(parents=True, exist_ok=True)
    py2 = conda_env_python(s.interpreter) if s.interpreter else None
    if s.interpreter and not py2:
        py2 = _create_companion_env(s.interpreter, s.interpreter_spec)
        if not py2:
            raise LncKitError(
                f"{tool}: needs the '{s.interpreter}' env (not found and could not be "
                f"created). Create it beside the lnckit env, then re-install {tool}.")

    exposed = []
    for want in s.exposes:
        found = _find_in_tree(srcroot, want)
        if not found:
            raise LncKitError(
                f"{tool}: '{want}' not found after install in {srcroot}")
        found.chmod(found.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        link = bindir / want
        if link.exists() or link.is_symlink():
            link.unlink()
        if py2:
            # launcher that resolves the interpreter env's python AT RUN TIME, so a
            # rebuilt/moved env (or a different envs dir than at install) does not leave
            # a dead absolute path behind. See _py2_shim_text.
            link.write_text(_py2_shim_text(s.interpreter, found, py2))
            link.chmod(0o755)
        else:
            link.symlink_to(found)
        exposed.append(want)
    for alias, target in (s.aliases or {}).items():
        link = bindir / alias
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(bindir / target)
        exposed.append(alias)

    ensure_managed_bin_on_path()
    return f"source: installed {', '.join(exposed)} -> {bindir}"


def build_conda_command(pkgs: list, env: Optional[str], yes: bool,
                        allow_missing: bool = False) -> list:
    exe = conda_exe()
    if not exe:
        if not allow_missing:
            raise LncKitError(
                "conda/mamba not found on PATH. Install Miniforge/Miniconda first, "
                "or use --method direct for tools that ship a prebuilt binary.")
        exe = "conda"  # placeholder for dry-run display
    name = os.path.basename(exe)
    # conda-forge MUST come before bioconda (bioconda packages depend on conda-forge),
    # and --override-channels drops the `defaults` channel - otherwise, with
    # `channel_priority: strict`, the solve fails or picks ancient builds.
    cmd = [name, "install", "--override-channels", "-c", "conda-forge", "-c", "bioconda"]
    if env:
        cmd += ["-n", env]
    if yes:
        cmd.append("-y")
    cmd += sorted(set(pkgs))
    return cmd


def install_via_conda(pkgs: list, env: Optional[str], yes: bool,
                      dry_run: bool = False) -> str:
    cmd = build_conda_command(pkgs, env, yes, allow_missing=dry_run)
    if dry_run:
        return "conda: would run: " + " ".join(cmd)
    import subprocess
    log.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)  # noqa: S603
    return "conda: " + " ".join(cmd)


# System package managers, in the order we prefer them, with their install verb.
_PKG_MANAGERS = [
    ("dnf", ["install", "-y"]),
    ("apt-get", ["install", "-y"]),
    ("yum", ["install", "-y"]),
    ("zypper", ["--non-interactive", "install"]),
    ("pacman", ["-S", "--noconfirm"]),
    ("brew", ["install"]),          # brew must not be run under sudo
]

# Build-time dependencies a source compile needs, per tool, as system package
# names (used to advise the user; installed via the same system package manager).
SYSTEM_BUILD_DEPS: dict[str, dict] = {
    "samtools": {
        "dnf":     ("gcc", "make", "zlib-devel", "bzip2-devel", "xz-devel",
                    "ncurses-devel", "libcurl-devel"),
        "apt-get": ("gcc", "make", "zlib1g-dev", "libbz2-dev", "liblzma-dev",
                    "libncurses5-dev", "libcurl4-openssl-dev"),
    },
}


def detect_pkg_manager() -> Optional[tuple]:
    """Return (exe, install_args) for the first available system package manager."""
    for exe, args in _PKG_MANAGERS:
        if shutil.which(exe):
            return exe, args
    return None


def _system_install(pkgs: list, dry_run: bool = False) -> str:
    """Install OS packages via the detected package manager (sudo for most)."""
    import subprocess
    mgr = detect_pkg_manager()
    if not mgr:
        raise LncKitError(
            "No system package manager found (dnf/apt-get/yum/zypper/pacman/brew).")
    exe, verb = mgr
    sudo = [] if exe == "brew" or os.geteuid() == 0 else ["sudo"]
    cmd = sudo + [exe, *verb, *pkgs]
    if dry_run:
        return "system: would run: " + " ".join(cmd)
    log.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)  # noqa: S603 - curated package names
    return "system: " + " ".join(cmd)


def install_via_system(pkg: str, dry_run: bool = False) -> str:
    """Install a single OS package (the `system` install method)."""
    return _system_install([pkg], dry_run)


def install_build_deps(tool: str, dry_run: bool = False) -> Optional[str]:
    """Install a source tool's compile dependencies via the system package manager."""
    mgr = detect_pkg_manager()
    if not mgr:
        return None
    deps = SYSTEM_BUILD_DEPS.get(tool, {}).get(mgr[0])
    if not deps:
        return None
    return _system_install(list(deps), dry_run)


def install_via_pip(pkg: str, dry_run: bool = False, break_system: bool = True) -> str:
    cmd = [sys.executable, "-m", "pip", "install", pkg]
    if break_system:
        cmd.append("--break-system-packages")
    if dry_run:
        return "pip: would run: " + " ".join(cmd)
    import subprocess
    log.info("running: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)  # noqa: S603
    except subprocess.CalledProcessError:
        # retry without the flag for normally-managed environments
        cmd2 = [sys.executable, "-m", "pip", "install", pkg]
        subprocess.run(cmd2, check=True)  # noqa: S603
        return "pip: " + " ".join(cmd2)
    return "pip: " + " ".join(cmd)


# --------------------------------------------------------------------------- #
# High-level planning + execution
# --------------------------------------------------------------------------- #
def choose_method(tool: str, prefer: str = "auto") -> Optional[str]:
    """Pick an install method for a tool given availability and preference."""
    recipe = INSTALL_RECIPES.get(tool)
    if not recipe:
        return None
    avail = recipe.methods()
    if not avail:
        return None
    if prefer != "auto" and prefer in avail:
        # 'direct' is only usable if we actually have a recipe for this platform
        return prefer
    have_conda = conda_exe() is not None
    # Preference order: conda (if usable) > pip > source > direct > system
    if "conda" in avail and have_conda:
        return "conda"
    if "pip" in avail:
        return "pip"
    if "source" in avail:
        return "source"
    if "direct" in avail:
        return "direct"
    if "system" in avail:
        return "system"
    if "conda" in avail:
        return "conda"  # will error with a helpful "install conda" message
    return avail[0]


@dataclass
class InstallPlan:
    conda_pkgs: list = field(default_factory=list)
    pip_pkgs: list = field(default_factory=list)
    direct_tools: list = field(default_factory=list)
    source_tools: list = field(default_factory=list)
    system_tools: list = field(default_factory=list)  # (tool, os_package) pairs
    unresolved: list = field(default_factory=list)   # tools with no method


def plan_install(missing: list, prefer: str = "auto") -> InstallPlan:
    plan = InstallPlan()
    for tool in missing:
        method = choose_method(tool, prefer)
        recipe = INSTALL_RECIPES.get(tool)
        if method == "conda":
            plan.conda_pkgs.append(recipe.conda)
        elif method == "pip":
            plan.pip_pkgs.append(recipe.pip)
        elif method == "source":
            plan.source_tools.append(tool)
        elif method == "direct":
            plan.direct_tools.append(tool)
        elif method == "system":
            plan.system_tools.append((tool, recipe.system))
        else:
            plan.unresolved.append(tool)
    return plan


def execute_plan(plan: InstallPlan, env: Optional[str] = None, yes: bool = False,
                 dry_run: bool = False) -> list:
    """Run a plan; returns a list of human-readable result lines."""
    results = []
    if plan.conda_pkgs:
        results.append(install_via_conda(plan.conda_pkgs, env, yes, dry_run))
    for pkg in plan.pip_pkgs:
        results.append(install_via_pip(pkg, dry_run))
    seen_sources = set()            # one shared archive (e.g. sra-tools) installs once
    for tool in plan.source_tools:
        recipe = INSTALL_RECIPES.get(tool)
        sig = (recipe.source.url or recipe.source.git) if recipe and recipe.source else tool
        if sig in seen_sources:
            results.append(f"source: {tool} already provided by an earlier install")
            continue
        seen_sources.add(sig)
        try:
            results.append(install_from_source(tool, dry_run))
        except Exception as e:      # noqa: BLE001 - clone/build can fail many ways
            results.append(f"source: FAILED for {tool}: {e}")
    for tool in plan.direct_tools:
        try:
            results.append(install_direct(tool, dry_run))
        except LncKitError as e:
            results.append(f"direct: FAILED for {tool}: {e}")
    for tool, pkg in plan.system_tools:
        try:
            results.append(install_via_system(pkg, dry_run))
        except Exception as e:      # noqa: BLE001 - sudo/pkg-mgr can fail many ways
            results.append(f"system: FAILED for {tool}: {e}")
    for tool in plan.unresolved:
        results.append(f"skip: no install recipe available for {tool}")
    return results
