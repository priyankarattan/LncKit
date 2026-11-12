#!/usr/bin/env python3
"""LncKit - one command to run the whole lncRNA identification pipeline.

The pipeline is launched with `lnckit.sh` (it sets up the tool environment, then runs
this script). `./lnckit.py` takes exactly the same arguments and works too. Driven
entirely by the command line - there is no interactive mode; set only the steps you
want to change and anything you leave out uses its default:
    ./lnckit.sh <folder>                   # run the default pipeline
    ./lnckit.sh <folder> --aligner star    # change one step; the rest stay default
    ./lnckit.sh <folder> --organism 1 --coding 1,2,5
Run  ./lnckit.sh --help  to see every option.
"""
import os
import sys
from pathlib import Path


def _bootstrap_env():
    """Make `./lnckit.py ...` work on its own by putting the pipeline's tools on PATH.

    The bioinformatics tools live in conda environments; lnckit.sh knows how to find them
    and set PATH/PYTHONPATH/etc. Rather than duplicate (and risk drifting from) that
    proven logic, when lnckit.py is started outside the environment we simply re-launch
    through lnckit.sh, which sets everything up and calls this script back. `lnckit.sh` exports
    `_LNCKIT_BOOTSTRAPPED=1` before doing so, which is how we avoid an infinite loop.
    """
    if os.environ.get("_LNCKIT_BOOTSTRAPPED"):
        return                                   # already inside the env - carry on
    launcher = Path(__file__).resolve().parent / "lnckit.sh"
    if launcher.exists():
        os.environ["_LNCKIT_BOOTSTRAPPED"] = "1"
        os.execv("/bin/bash", ["/bin/bash", str(launcher), *sys.argv[1:]])
    # No lnckit.sh next to us: assume the environment is already active (the imports below
    # will raise a clear error if it is not) and continue.


if __name__ == "__main__":
    _bootstrap_env()

import argparse   # noqa: E402  (imported after the bootstrap above, on purpose)
import csv        # noqa: E402
import re         # noqa: E402

from lnckit.core.config import CHOICES, SPLICE_AWARE, load_config  # noqa: E402
from lnckit.core.exceptions import ConfigError                     # noqa: E402

# The steps the user chooses tools for, in pipeline order.
STAGES = [
    ("qc",               "Quality control"),
    ("trim",             "Trimming"),
    ("aligner",          "Alignment"),
    ("assembler",        "Assembly"),
    ("coding_potential", "Coding potential"),
    ("protein_filter",   "Protein homology"),
    ("housekeeping",     "ncRNA filtering"),
    ("classifier",       "Classification"),
]


# Non-interactive (batch) mode: no terminal attached (e.g. run under qsub/sbatch/nohup)
# or forced with LNCKIT_BATCH=1. In this mode every prompt takes its default WITHOUT
# reading stdin, so the run is deterministic and safe for schedulers - configure the
# whole run via samples.txt + pipeline.yaml instead of answering questions.
BATCH = (os.environ.get("LNCKIT_BATCH") == "1") or (not sys.stdin.isatty())


def prompt(text, default=""):
    """input() that survives being piped/EOF and falls back to a default.

    In BATCH mode we never touch stdin: we echo the prompt with the default it took,
    so the job log shows exactly what each step chose.
    """
    if BATCH:
        print(f"{text}{default}   [batch: default]")
        return default
    try:
        return input(text).strip() or default
    except EOFError:
        return default


_ORGANISMS = [
    ("plant",     "Plant      (CPAT plant model, CNCI 'pl')"),
    ("human",     "Human      (CPAT cutoff 0.364, CNCI 've')"),
    ("mouse",     "Mouse      (CPAT cutoff 0.44,  CNCI 've')"),
    ("fly",       "Fly / Drosophila"),
    ("zebrafish", "Zebrafish"),
    ("other",     "Other / unknown  (CPAT is SKIPPED - no model applied)"),
]


_STRANDNESS = [
    ("auto", "Auto-detect from the reads        - default; falls back to unstranded if unsure"),
    ("RF",   "Stranded, dUTP first-strand (RF)  - most common; StringTie --rf"),
    ("FR",   "Stranded, second-strand     (FR)  - StringTie --fr"),
    ("none", "Unstranded                        - no --rf/--fr (only if truly unstranded)"),
]

# Which command-line flag drives each pipeline stage (some read nicer than the config key).
_STEP_FLAG = {
    "qc": "qc", "trim": "trim", "aligner": "aligner", "assembler": "assembler",
    "coding_potential": "coding", "protein_filter": "protein",
    "housekeeping": "housekeeping", "classifier": "classifier",
}


def _url_target(url, key):
    """Local filename to save a reference URL to (basename, minus any .gz/query)."""
    name = url.rstrip("/").split("/")[-1].split("?")[0] or key
    return name[:-3] if name.endswith(".gz") else name


def _store_reference(ref, ref_key, val):
    """Record one reference value from the sheet.

    Accepts a local PATH or a URL, interchangeably, for any reference. A URL (http/https/
    ftp, or an explicit `<key>_url` row) is stored as `<key>_url` plus a derived local
    target path, so setup_reference() downloads it; a plain value is stored as the path.
    """
    if ref_key.endswith("_url"):
        ref[ref_key] = val
        ref.setdefault(ref_key[:-4], _url_target(val, ref_key[:-4]))
    elif val.lower().startswith(("http://", "https://", "ftp://")):
        ref[ref_key + "_url"] = val
        ref.setdefault(ref_key, _url_target(val, ref_key))
    else:
        ref[ref_key] = val


_REF_KEYS = {
    "genome_fasta": "genome_fasta", "genome": "genome_fasta", "fasta": "genome_fasta",
    "annotation_gtf": "annotation_gtf", "gtf": "annotation_gtf",
    "annotation_gff": "annotation_gff", "gff": "annotation_gff",
    "annotation_gff3": "annotation_gff", "gff3": "annotation_gff",
    "protein_db": "protein_db", "rfam_cm": "rfam_cm",
    "genome_protein": "genome_protein", "species_protein": "genome_protein",
    "proteins": "genome_protein", "pep": "genome_protein",
}


#: sheet keys that name a LOCAL folder of FASTQ files (turns the run into a local run).
_READS_DIR_KEYS = {"reads_dir", "reads", "reads_folder", "fastq_dir", "fastq"}


def read_project_sheet(path="samples.txt"):
    """Parse the one-file project sheet: genome/reference paths + samples (+ local reads).

    A single text sheet (comma-separated) describes the whole run. In any order it may contain:
      * reference rows -> `<key>,<path>` where key is genome_fasta (or genome/fasta),
        annotation_gtf (or gtf), annotation_gff (or gff), protein_db, rfam_cm.
      * a local-reads row -> `reads_dir,<folder>` (makes this a LOCAL run: reads are read
        from that folder instead of being downloaded from SRA).
      * sample rows -> `<name>,<group>`, optionally preceded by an `accession,group`
        (or `sample,group`) header line. For an SRA run <name> is a run accession; for a
        local run <name> is the FASTQ file prefix.
    Blank lines and `#` comments are ignored.

    Returns (accessions, {group: [names]}, reference_overrides, reads_dir_or_None).
    """
    accs, groups, ref, reads_dir = [], {}, {}, None
    if not Path(path).exists():
        return accs, groups, ref, reads_dir
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            first = row[0].strip()
            if not first or first.startswith("#"):
                continue
            key = first.lower()
            if key in ("accession", "accessions", "sample", "srr"):
                continue                                  # sample-section header
            if key in _READS_DIR_KEYS:                    # local reads folder
                if len(row) > 1 and row[1].strip():
                    reads_dir = row[1].strip()
                continue
            # A reference row: `<key>,<path-or-URL>`. A ref key is NEVER a sample row -
            # even with a blank value - so e.g. `genome_protein` can't fall through and
            # get fetched as an SRA accession (which 404'd in prefetch/fastq-dump).
            ref_key = _REF_KEYS.get(key)
            if ref_key is None and key.endswith("_url") and key[:-4] in _REF_KEYS:
                ref_key = _REF_KEYS[key[:-4]] + "_url"    # explicit `<key>_url` row
            if ref_key is not None:
                val = row[1].strip() if len(row) > 1 else ""
                if val:
                    _store_reference(ref, ref_key, val)
                continue
            acc = first                                   # otherwise: a sample row
            grp = row[1].strip() if len(row) > 1 and row[1].strip() else ""
            accs.append(acc)
            groups.setdefault(grp, []).append(acc)
    return accs, groups, ref, reads_dir


def run_config_path(output_dir):
    """Where the resolved config of the current run is saved for later resume."""
    return Path(output_dir) / ".lnckit" / "run_config.yaml"


def save_run_config(cfg):
    """Persist the resolved config so a later run can resume without re-asking."""
    path = run_config_path(cfg.output_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(cfg.to_yaml())


def has_progress(output_dir):
    """True if a previous run left any resumable output behind (reads/BAMs/assembly)."""
    out = Path(output_dir)
    for sub in ("raw_reads", "trimmed_reads", "alignment", "assembly"):
        d = out / sub
        if d.is_dir() and any(d.iterdir()):
            return True
    return False


def setup_tools(cfg):
    """Check every tool this pipeline needs; install the missing ones.

    Prints one line per tool: [OK] if present, [MISSING] otherwise, then installs
    the missing ones and confirms. Returns True if all tools are present at the end.
    """
    from lnckit.core.installer import doctor, execute_plan, plan_install

    print("\n" + "-" * 52)
    print("STEP 1 - required tools")
    print("-" * 52)
    statuses = doctor(cfg)
    missing = []
    for s in statuses:
        if s.present:
            print(f"   [OK]      {s.tool}" + (f"  (v{s.version})" if s.version else ""))
        else:
            print(f"   [MISSING] {s.tool}")
            missing.append(s.tool)

    if not missing:
        print("\n   All required tools are already installed.")
        return True

    print(f"\n   {len(missing)} tool(s) to install: {', '.join(missing)}")
    if prompt("   Install them now? [Y/n]: ", "y").lower() == "n":
        return False
    print("   Installing (this can take a few minutes)...")
    # Install conda packages INTO the env this interpreter runs from (the lnckit env,
    # which is what lnckit.sh puts on PATH) - not the base env, which lnckit.sh never sees.
    env_name = os.path.basename(sys.prefix) or None
    for line in execute_plan(plan_install(missing), env=env_name, yes=True):
        print("     ", line)

    still = [s.tool for s in doctor(cfg) if not s.present]
    if still:
        print(f"\n   Still missing: {', '.join(still)}")
        # CNCI/PLEK are NOT conda packages - they are Python-2 source tools that need
        # the lnckit_py2 env. On a shared/read-only conda base, create it in your own
        # env dir, then re-run (LncKit will build CNCI/PLEK against it).
        if {"cnci", "plek"} & set(still):
            print("   CNCI/PLEK need a Python-2 env. Create it once, then re-run:")
            print("     conda create -n lnckit_py2 --override-channels "
                  "-c conda-forge -c bioconda -y python=2.7 biopython numpy")
        other = [t for t in still if t not in ("cnci", "plek")]
        if other:
            print("   For the rest, try:  conda install -n lnckit -c conda-forge "
                  "-c bioconda " + " ".join(other))
        return False
    print("\n   Dependencies completed - all tools present.")
    return True


def setup_reference(cfg):
    """Check the genome + annotation exist; download from a URL if one is configured.

    Add e.g. `genome_fasta_url: https://...` under `reference:` in pipeline.yaml to
    have it downloaded automatically (gzipped files are unpacked).
    """
    import gzip
    import shutil
    import urllib.request

    ref = cfg.get("reference", {}) or {}
    print("\n" + "-" * 52)
    print("STEP 2 - reference files (genome / GTF)")
    print("-" * 52)
    ok = True
    for key in ("genome_fasta", "annotation_gtf", "annotation_gff", "protein_db",
                "genome_protein", "rfam_cm", "cpat_hexamer", "cpat_logit"):
        path = ref.get(key)
        if isinstance(path, (list, tuple)):          # protein_db can be a list of DBs
            for i, p in enumerate(path):
                mark = "[OK]     " if (p and Path(p).exists()) else "[MISSING]"
                print(f"   {mark} {key}[{i}]: {p}")
            continue
        needed = key in ("genome_fasta", "annotation_gtf")
        if not path:
            # An unset path is fine for the optional files; for the genome/annotation it
            # means none was configured (no path or URL given).
            if needed:
                print(f"   [MISSING] {key}: not set - give a path or URL in samples.txt "
                      f"(or pipeline.yaml under reference:).")
                ok = False
            continue
        if Path(path).exists():
            print(f"   [OK]      {key}: {path}")
            continue
        url = ref.get(f"{key}_url")
        if url:
            print(f"   [DOWNLOAD] {key}: fetching...")
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            # Download + unpack to a temp file, then move into place only on success,
            # so an interrupted run never leaves a partial file masquerading as done.
            dl = str(path) + (".gz.part" if url.endswith(".gz") else ".part")
            part = str(path) + ".part"
            try:
                urllib.request.urlretrieve(url, dl)
                if url.endswith(".gz"):
                    with gzip.open(dl, "rb") as fi, open(part, "wb") as fo:
                        shutil.copyfileobj(fi, fo)
                    Path(dl).unlink()
                os.replace(part, path)          # atomic on the same filesystem
                print(f"   [OK]      {key}: -> {path}")
            except BaseException:
                for leftover in (dl, part):     # clean partials on error/interrupt
                    Path(leftover).unlink(missing_ok=True)
                raise
        else:
            print(f"   [{'MISSING' if needed else 'skip   '}] {key}: put file at {path}"
                  f"  (or add {key}_url in pipeline.yaml)")
            if needed:
                ok = False
    # Homology DBs auto-included at the protein step if present in $LNCKIT_DB_DIR (COG /
    # Pfam-A downloaded via download_databases.sh). Shown here so it's clear they'll be used
    # even though they aren't reference.* config keys.
    _dbdir = Path(os.environ.get("LNCKIT_DB_DIR", "databases"))
    for _lbl, _fn in (("cog-20.fa", "cog-20.fa"), ("Pfam-A.fasta", "Pfam-A.fasta")):
        _f = _dbdir / _fn
        if _f.exists():
            print(f"   [OK]      protein_db (auto): {_f}")
    if not ok:
        print("\n   >> How to provide your genome files (ANY filename is fine):")
        print( "      List them in samples.txt as reference rows - a local PATH or a URL:")
        print( "          genome_fasta,/path/to/YourGenome.fasta")
        print( "          annotation_gtf,/path/to/YourGenome.gtf     # a .gff / .gff3 also works")
        print( "      A URL is auto-downloaded, e.g.:")
        print( "          genome_fasta,https://host/YourGenome.fasta.gz")
        print( "      Or set reference.genome_fasta / annotation_gtf in pipeline.yaml.")
    return ok


def _ref_summary(ref):
    return ("; reference from sheet: "
            + ", ".join(f"{k}={v}" for k, v in ref.items())) if ref else ""


def _step_range():
    """Resolve --from/--until (passed by lnckit.sh via env) to (from_step, until_step).

    Lets a run JUMP to a step - e.g. you already have SAM/BAM files:
        lnckit.sh <dir> --from align       # skip fetch/qc/trim/index, start at alignment
    Validated against the real step list so a typo fails fast, not mid-run.
    """
    from lnckit.workflow.steps import STEP_NAMES
    fs = os.environ.get("LNCKIT_FROM_STEP") or None
    us = os.environ.get("LNCKIT_UNTIL_STEP") or None
    for label, v in (("--from", fs), ("--until", us)):
        if v and v not in STEP_NAMES:
            print(f"Unknown step {v!r} for {label}. Valid steps:\n  {', '.join(STEP_NAMES)}")
            sys.exit(2)
    if fs or us:
        print(f"Step range: from={fs or 'start'}  until={us or 'end'}")
    return fs, us


def _launch_engine(cfg):
    """Run the pipeline, turning a tool failure into a concise message (not a traceback)."""
    from lnckit.workflow.engine import Engine
    from lnckit.core.exceptions import ToolExecutionError, OutputValidationError
    fs, us = _step_range()
    try:
        Engine(cfg).run(from_step=fs, until_step=us)
    except (ToolExecutionError, OutputValidationError) as e:
        # A tool failed mid-run. The step and command are already in the log; show a
        # concise message (not a Python traceback) and stop with a non-zero status.
        print("\n" + "=" * 52, file=sys.stderr)
        print("RUN STOPPED - a pipeline step failed:", file=sys.stderr)
        print(f"   {e}", file=sys.stderr)
        print(f"\nFull log: {Path(cfg.output_dir) / 'logs' / 'lnckit.log'}", file=sys.stderr)
        print("Fix the cause above, then re-run the same folder to resume "
              "from where it stopped.", file=sys.stderr)
        sys.exit(1)


def run_resume(cfg):
    """Skip straight to running an existing run's saved config, continuing from the
    step it stopped at. Tools/reference are still verified (fast when present); the
    engine reuses every finished step's outputs and re-does only what is missing."""
    print("\n" + "-" * 52)
    print("RESUME - continuing the previous run where it stopped")
    print("-" * 52)
    accs = cfg.get("inputs.sra_accessions", []) or []
    print(f"Samples : {len(accs)} accession(s)")
    print("Tools   :")
    for key, name in STAGES:
        val = cfg.tools.get(key)
        if isinstance(val, list):
            if key == "coding_potential":
                eff = cfg.get("coding_potential.consensus", "all")
                val = " + ".join(val) + f"   ({_consensus_label(eff)})"
            else:
                val = " + ".join(val) + "   (transcript must pass ALL)"
        print(f"   {name:18s}: {val}")

    if not setup_tools(cfg):
        print("\nStopping: not all required tools are available.")
        return
    if not setup_reference(cfg):
        print("\nStopping: genome / annotation not available (see above).")
        return
    print("\n" + "-" * 52)
    print("Resuming the pipeline...")
    print("-" * 52)
    _launch_engine(cfg)


# =====================================================================
# Command-line interface
# ---------------------------------------------------------------------
# LncKit is driven entirely by the command line - there is no interactive mode.
# Every step is a flag; set only what you want to change and the rest use their
# defaults, so the whole run can be launched from a script and left on an HPC node
# (nohup / sbatch / qsub) with no terminal open. Each step flag takes a 1-based
# NUMBER or the tool's name. With no flags at all, the default pipeline runs.
# =====================================================================
_STEP_FLAG_KEYS = ["qc", "trim", "aligner", "assembler", "coding",
                   "protein", "housekeeping", "classifier"]


def _consensus_label(mode):
    """Human description of how a coding-potential ensemble combines its tools' votes."""
    m = str(mode or "all").strip().lower()
    if m in ("majority", "consensus", "vote"):
        return "kept if a majority call it non-coding"
    if m.isdigit():
        return f"kept if >= {m} tools call it non-coding"
    return "kept only if ALL tools call it non-coding"


def _resolve_menu(value, pairs, flag):
    """Map a CLI value (a 1-based menu number OR a name) to its canonical name."""
    names = [p[0] if isinstance(p, tuple) else p for p in pairs]
    v = str(value).strip()
    if v.isdigit():
        n = int(v)
        if not (1 <= n <= len(names)):
            raise SystemExit(f"ERROR: --{flag} {n} is out of range (choose 1..{len(names)}).")
        return names[n - 1]
    lower = [n.lower() for n in names]                 # match names case-insensitively
    if v.lower() in lower:
        return names[lower.index(v.lower())]           # ...but return the canonical spelling
    raise SystemExit(f"ERROR: --{flag} '{value}' invalid - use 1..{len(names)} or one of {names}.")


_USAGE = (
    "lnckit.sh [-h] [--input PATH] [--paired | --single] [--strandness N]\n"
    "                 [--merge | --no-merge] [--qc N] [--trim N] [--aligner N]\n"
    "                 [--assembler N] [--organism N] [--coding LIST] [--consensus V]\n"
    "                 [--protein N] [--housekeeping N] [--classifier N]\n"
    "                 [--output DIR] [--config FILE] [--threads N]\n"
    "                 [--from STEP] [--until STEP] [--dry-run] [run-folder]")


def _cli_help_text():
    """The full, ordered `--help` body: usage guide + one table with a DEFAULT column."""
    d = load_config(None)                        # packaged defaults, for the DEFAULT column

    def opts(pairs):
        return " ".join(f"{i}={p[0] if isinstance(p, tuple) else p}"
                        for i, p in enumerate(pairs, 1))

    def row(flag, choices, default=""):
        return f"   {flag:<22}{choices:<52}{default}"

    strand_def = d.get("inputs.strandness") or "auto"
    L = [
        "=" * 64,
        "   LncKit - lncRNA identification pipeline",
        "=" * 64,
        "",
        "USAGE",
        "   ./lnckit.sh <run-folder> [options]",
        "",
        "   Run the default pipeline .....  ./lnckit.sh .",
        "   Change only what you need ....  ./lnckit.sh . --aligner star --organism plant",
        "   Use local reads (not SRA) ....  ./lnckit.sh . --input reads/",
        "   Preview the plan, don't run ..  ./lnckit.sh . --dry-run",
        "",
        "   <run-folder> holds samples.txt (+ genome files); results go to",
        "   <run-folder>/lnckit_output/. Any step you don't set uses its DEFAULT.",
        "   Give a tool by NUMBER or NAME.  (./lnckit.py takes the same arguments.)",
        "",
        "INPUT",
        row("FLAG", "CHOICES / VALUE", "DEFAULT"),
        row("-i, --input PATH", "samples.txt (sheet)  OR  a folder of reads", "samples.txt"),
        row("--paired / --single", "read layout (otherwise auto-detected)", "auto"),
        row("-str, --strandness", opts(_STRANDNESS), strand_def),
        row("--merge / --no-merge", "pool replicate runs by group", "merge (on)"),
        "",
        "PIPELINE STEPS  (in order; give a tool by number or name)",
        row("--qc", opts(CHOICES["qc"]), d.tools.get("qc")),
        row("--trim", opts(CHOICES["trim"]), d.tools.get("trim")),
        row("-a, --aligner", opts(CHOICES["aligner"]), d.tools.get("aligner")),
        row("--assembler", opts(CHOICES["assembler"]), d.tools.get("assembler")),
        row("-org, --organism", opts(_ORGANISMS), d.get("organism", "other")),
        row("--coding", opts(CHOICES["coding_potential"]), "1,2,3,4,5,6"),
        row("", "(comma-list several tools = an ensemble)"),
        row("--consensus", "all | majority | N",
            d.get("coding_potential.consensus", "majority")),
        row("--protein", opts(CHOICES["protein_filter"]), d.tools.get("protein_filter")),
        row("--housekeeping", opts(CHOICES["housekeeping"]), d.tools.get("housekeeping")),
        row("--classifier", opts(CHOICES["classifier"]), d.tools.get("classifier")),
        "",
        "RUN CONTROL",
        row("-o, --output DIR", "where results go", "<run-folder>/lnckit_output"),
        row("-c, --config FILE", "advanced YAML settings", "pipeline.yaml (if present)"),
        row("-p, --threads N", "CPU threads", "8"),
        row("--from / --until STEP", "run only part of the pipeline"),
        row("--dry-run", "print the plan and exit without running"),
    ]
    return "\n".join(L)


def build_cli_parser():
    # The whole help body lives in `description` (printed verbatim); every option's help is
    # suppressed so argparse does not append a second, duplicated option list. A custom
    # `usage=` keeps the compact synopsis line at the very top.
    p = argparse.ArgumentParser(
        prog="lnckit.sh", usage=_USAGE, description=_cli_help_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=True)
    S = argparse.SUPPRESS
    p.add_argument("run_folder", nargs="?", help=S)   # a directory (results go there); a file = sheet
    p.add_argument("-i", "--input", metavar="PATH", help=S)  # samples.txt (sheet) OR a reads folder
    p.add_argument("-org", "--organism", dest="organism", help=S)
    p.add_argument("-str", "--strandness", "--strand", dest="strandness", help=S)
    p.add_argument("--qc", help=S)
    p.add_argument("--trim", help=S)
    p.add_argument("-a", "--aligner", "--align", dest="aligner", help=S)
    p.add_argument("--assembler", "--assembly", dest="assembler", help=S)
    p.add_argument("--coding", "--coding-potential", dest="coding", help=S)
    p.add_argument("--protein", "--protein-filter", dest="protein", help=S)
    p.add_argument("--housekeeping", help=S)
    p.add_argument("--classifier", help=S)
    p.add_argument("--consensus", help=S)
    p.add_argument("--paired", action="store_true", help=S)
    p.add_argument("--single", action="store_true", help=S)
    p.add_argument("--merge", action="store_true", help=S)
    p.add_argument("--no-merge", dest="no_merge", action="store_true", help=S)
    p.add_argument("-o", "--output", help=S)
    p.add_argument("-c", "--config", help=S)
    p.add_argument("--from", dest="from_step", help=S)
    p.add_argument("--until", dest="until_step", help=S)
    p.add_argument("-p", "--threads", type=int, help=S)
    p.add_argument("--dry-run", dest="dry_run", action="store_true", help=S)
    return p


def _flags_given(args):
    """True if the user supplied any step / input flag (drives the intro message only)."""
    if args.dry_run:
        return True
    if any(getattr(args, f) for f in _STEP_FLAG_KEYS):
        return True
    if args.organism or args.strandness or args.input or args.consensus:
        return True
    return bool(args.paired or args.single or args.merge or args.no_merge)


def _merge_flag(args, groups):
    """Pool replicate runs by group. ON by default; --no-merge disables, --merge forces."""
    if args.no_merge:
        return False
    if args.merge:
        return True
    # default (on): pool whenever a group holds more than one run.
    return bool(groups) and any(len(a) > 1 for a in groups.values())


def _detect_layout(reads_dir, args):
    """Decide paired- vs single-end for a LOCAL reads folder.

    Auto-detects from the filenames (``*_1.fastq*`` with ``*_2.fastq*`` -> paired), but an
    explicit --paired / --single always wins. Returns (paired, n_samples, note_string).
    """
    rd = Path(reads_dir)
    r1 = sorted(rd.glob("*_1.f*q*"))
    singles = [p for p in rd.glob("*.f*q*")
               if "_1." not in p.name and "_2." not in p.name]
    if args.single:
        paired, how = False, "set by --single"
    elif args.paired:
        paired, how = True, "set by --paired"
    elif r1:
        paired, how = True, "auto-detected: found *_1 / *_2 pairs"
    elif singles:
        paired, how = False, "auto-detected: no paired *_1/*_2 files"
    else:
        paired, how = True, "no FASTQ matched"
    n = len(r1) if paired else len(singles)
    note = (f"Reads    : {n} {'paired-end' if paired else 'single-end'} "
            f"sample(s) in {rd}   [{how}]")
    return paired, n, note


def _tools_from_flags(cfg_file, args):
    """Start from the config defaults and override only the stages named on the CLI."""
    tools = dict(load_config(cfg_file).tools)
    for key, _name in STAGES:
        val = getattr(args, _STEP_FLAG.get(key, key), None)
        if val is None:
            continue
        opts = CHOICES[key]
        if key == "coding_potential":
            picks = [_resolve_menu(x, opts, "coding")
                     for x in re.split(r"[,\s]+", val.strip()) if x]
            if picks:
                tools[key] = picks if len(picks) > 1 else picks[0]
        else:
            tools[key] = _resolve_menu(val, opts, _STEP_FLAG.get(key, key))
    return tools


def _input_from_flags(args):
    """Build the inputs override (SRA download OR local FASTQ) from --input + the sheet.

    --input may be a project sheet (a file) or a folder of local reads (a directory); if
    omitted it defaults to samples.txt. A sheet carrying a `reads_dir` row is a local run
    too. Local runs auto-detect paired/single unless --paired/--single is given.
    """
    inp = args.input
    sheet, reads_dir = "samples.txt", None
    if inp:
        if Path(inp).expanduser().is_dir():
            reads_dir = inp                       # --input points at a reads folder
        else:
            sheet = inp                           # --input points at a project sheet
    accessions, groups, ref, sheet_reads = (
        read_project_sheet(sheet) if Path(sheet).exists() else ([], {}, {}, None))
    reads_dir = reads_dir or sheet_reads          # --input folder wins over a sheet reads_dir

    if reads_dir:                                 # ---------- LOCAL FASTQ run ----------
        rd = Path(reads_dir).expanduser()
        if not rd.is_dir():
            raise SystemExit(f"ERROR: reads folder not found: {reads_dir}")
        paired, n, note = _detect_layout(rd, args)
        print(note)
        if n == 0:
            raise SystemExit(
                f"ERROR: no FASTQ files in {rd} (expected <sample>_1.fastq[.gz] + "
                "<sample>_2.fastq[.gz] for paired-end, or <sample>.fastq[.gz] for "
                "single-end; .fastq, .fastq.gz, .fq and .fq.gz are all accepted).")
        merge = _merge_flag(args, groups)
        summary = (f"{n} local {'paired-end' if paired else 'single-end'} sample(s) from {rd}"
                   + ("  [pooled by group]" if merge else "") + _ref_summary(ref))
        overrides = {"inputs": {"reads_dir": str(rd), "sra_accessions": None,
                                "paired": paired, "groups": groups,
                                "merge_replicates": merge}}
    else:                                         # ---------- SRA download run ----------
        if not accessions:
            raise SystemExit(
                f"ERROR: no input found. Put SRA accessions in '{sheet}' (accession,group "
                "rows), or point --input at a folder of FASTQ files (or a sheet with a "
                "reads_dir row) for local data.")
        paired = not args.single                  # SRA: paired by default, --single overrides
        merge = _merge_flag(args, groups)
        label = ", ".join(f"{g} ({len(a)})" for g, a in groups.items()) or "(none)"
        summary = (f"{len(accessions)} {'paired-end' if paired else 'single-end'} "
                   f"SRA accession(s), groups: {label}"
                   + ("  [pooled by group]" if merge else "") + _ref_summary(ref))
        overrides = {"inputs": {"sra_accessions": accessions, "paired": paired,
                                "groups": groups, "merge_replicates": merge}}
    if ref:
        overrides["reference"] = ref
    return overrides, summary


def _reference_from_flags(cfg_file, src_overrides):
    """Genome/annotation come from samples.txt reference rows (a path or URL); pipeline.yaml
    is the fallback. No name-based guessing - the user always gives an explicit path or URL."""
    sheet_ref = (src_overrides or {}).get("reference", {})
    return dict(sheet_ref)


def _state_from_flags(cfg_file, args, env_out):
    """Build the run state (tools, input, reference, organism) from command-line flags."""
    tools = _tools_from_flags(cfg_file, args)
    src_overrides, input_summary = _input_from_flags(args)
    ref_over = _reference_from_flags(cfg_file, src_overrides)
    forced = (tools.get("aligner") not in SPLICE_AWARE
              and tools.get("assembler") in {"stringtie", "cufflinks"})
    organism = os.environ.get("LNCKIT_ORGANISM") \
        or load_config(cfg_file).get("organism", "other")
    return {
        "organism": organism,
        "tools": tools,
        "forced_unspliced": forced,
        "consensus_mode": args.consensus,          # None -> config default consensus rule
        "src_overrides": src_overrides,
        "input_summary": input_summary,
        "strandness": "__ASK__",                    # LNCKIT_STRANDNESS env is applied in load_config
        "ref_over": ref_over,
    }


def main():
    parser = build_cli_parser()
    # Bare `lnckit.py` (no arguments) shows the help instead of launching a run. Through
    # lnckit.sh this is handled there (it sees the real argument count, before the working
    # directory is consumed); this covers a direct `python lnckit.py` with no launcher.
    _bare = len(sys.argv) == 1 and not os.environ.get("_LNCKIT_BOOTSTRAPPED")
    _help = _bare or any(a in ("-h", "--help") for a in sys.argv[1:])
    if not _help:                                # skip the banner when only showing help
        print("=" * 52)
        print("   LncKit - lncRNA identification pipeline")
        print("=" * 52)
    if _bare:
        parser.print_help()
        return

    # ---- parse the command line. Every step is a flag; anything not set uses its
    #      default, and with no flags at all the default pipeline runs. ---------------
    args, _unknown = parser.parse_known_args()
    if _unknown:
        print(f"NOTE: ignoring unrecognised argument(s): {' '.join(_unknown)}")
    # Resolve path-like flags to ABSOLUTE before any chdir below, so a relative path on
    # the command line keeps pointing where it was typed (mirrors lnckit.sh). If lnckit.sh
    # already handled these (values in the env), leave them alone. (--input is intentionally
    # NOT resolved here: a relative --input should point inside the run-folder we cd into.)
    if args.config and not os.environ.get("LNCKIT_CONFIG"):
        args.config = str(Path(args.config).resolve())
    if args.output and not os.environ.get("LNCKIT_OUTPUT_DIR"):
        args.output = str(Path(args.output).resolve())
    # The positional is the working directory (like lnckit.sh's WORKDIR): a DIRECTORY is cd'd
    # into so results land there; a FILE is taken as the --input project sheet. (When the
    # call came through lnckit.sh the working dir was already consumed, so this is a no-op.)
    if args.run_folder:
        rf = Path(args.run_folder)
        if rf.is_dir():
            os.chdir(rf)
        elif rf.is_file() and not args.input:
            args.input = str(rf.resolve())
    # Flags that lnckit.sh also maps to env vars - set them here too so a DIRECT
    # `python lnckit.py ...` call behaves identically. lnckit.sh's values win if already set.
    for _name, _val in (("LNCKIT_CONFIG", args.config), ("LNCKIT_OUTPUT_DIR", args.output),
                        ("LNCKIT_ORGANISM", args.organism), ("LNCKIT_STRANDNESS", args.strandness),
                        ("LNCKIT_FROM_STEP", args.from_step), ("LNCKIT_UNTIL_STEP", args.until_step)):
        if _val is not None and not os.environ.get(_name):
            os.environ[_name] = str(_val)
    # Accept a menu NUMBER for organism/strandness wherever they arrive (flag or env),
    # resolving to the canonical name the rest of the code (and load_config) expects.
    if os.environ.get("LNCKIT_ORGANISM"):
        os.environ["LNCKIT_ORGANISM"] = _resolve_menu(
            os.environ["LNCKIT_ORGANISM"], _ORGANISMS, "organism")
    if os.environ.get("LNCKIT_STRANDNESS"):
        os.environ["LNCKIT_STRANDNESS"] = _resolve_menu(
            os.environ["LNCKIT_STRANDNESS"], _STRANDNESS, "strandness")

    # LncKit is driven entirely by the command line - there is no interactive mode.
    # Set only the steps you want; anything you leave out uses its default.
    if _flags_given(args):
        print("Running from command-line parameters. Any step you did not set uses its default.")
    else:
        print("No parameters given - running the DEFAULT pipeline. "
              "Run  ./lnckit.sh --help  to see every option.")

    # Config precedence: --config/-c (LNCKIT_CONFIG, an optional custom YAML) wins;
    # else a pipeline.yaml in the working dir; else the packaged defaults. Most users need
    # none of these - flags (or the defaults) drive the run.
    cfg_file = os.environ.get("LNCKIT_CONFIG") \
        or ("pipeline.yaml" if Path("pipeline.yaml").exists() else None)
    if cfg_file and not Path(cfg_file).exists():
        sys.exit(f"ERROR: config file not found: {cfg_file}")

    # ---- 0. resume a previous run? -----------------------------------
    # If an earlier run saved its config and left resumable outputs behind, offer to
    # pick up where it stopped instead of re-asking every question.
    # `LNCKIT_OUTPUT_DIR` (set by lnckit.sh -o/--output) overrides where results go.
    env_out = os.environ.get("LNCKIT_OUTPUT_DIR") or None
    output_dir = env_out or load_config(cfg_file).output_dir
    saved = run_config_path(output_dir)
    if saved.exists() and has_progress(output_dir):
        print(f"\nA previous run was found under '{output_dir}'.")
        if prompt("Resume it from where it stopped? [Y/n]: ", "y").lower() != "n":
            try:
                run_resume(load_config(saved))
            except ConfigError as e:
                print(f"\nSaved run config is unusable ({e}); starting fresh.")
            else:
                return

    # Build the whole run from the command line: every step you did not set with a flag
    # uses its default (from pipeline.yaml, else the packaged defaults). No questions.
    st = _state_from_flags(cfg_file, args, env_out)

    organism        = st["organism"]
    tools           = st["tools"]
    forced_unspliced = st["forced_unspliced"]
    consensus_mode  = st["consensus_mode"]
    src_overrides   = st["src_overrides"]
    input_summary   = st["input_summary"]
    strandness      = st["strandness"]
    ref_over        = st["ref_over"]

    # ---- assemble the config overrides from the answers collected above -----------
    overrides = {"tools": tools, **src_overrides}      # inputs + optional reference
    if organism:
        overrides["organism"] = organism
    if strandness != "__ASK__":                        # a value was chosen at the prompt
        overrides.setdefault("inputs", {})["strandness"] = strandness
    if env_out:
        overrides["output_dir"] = env_out              # -o/--output wins
    if args.threads:
        overrides.setdefault("resources", {})["threads"] = args.threads
    base_ref = dict((load_config(cfg_file).get("reference", {}) or {}))
    sheet_ref = overrides.get("reference", {})         # paths given in the sample sheet
    if ref_over:
        overrides["reference"] = ref_over

    def _ref_source(key):
        if sheet_ref.get(key):  return sheet_ref[key], "from samples.txt"
        if base_ref.get(key):   return base_ref[key],  "from pipeline.yaml"
        if ref_over.get(key):   return ref_over[key],  "you entered"
        return None, ""

    # ---- 5. show the plan - including EXACTLY which genome will be used ----------
    print("\n" + "-" * 52)
    print("PLAN")
    print("-" * 52)
    print(f"Input   : {input_summary}")
    print("Genome  :")
    for key, label in (("genome_fasta", "genome FASTA"), ("annotation_gtf", "annotation")):
        val, src = _ref_source(key)
        if val:
            print(f"   {label:12s}: {Path(val).resolve()}   ({src})")
        else:
            print(f"   {label:12s}: NOT SET - the run will stop until you provide it")
    if organism:
        print(f"Organism: {organism}   (CPAT/CNCI/LncFinder model)")
    print("Tools   :")
    for key, name in STAGES:
        val = tools.get(key)
        if isinstance(val, list):
            if key == "coding_potential":
                eff = consensus_mode or load_config(cfg_file).get(
                    "coding_potential.consensus", "all")
                val = " + ".join(val) + f"   ({_consensus_label(eff)})"
            else:
                val = " + ".join(val) + "   (transcript must pass ALL)"
        print(f"   {name:18s}: {val}")

    if forced_unspliced:
        # user accepted a non-splice-aware aligner; let validate() allow it
        overrides["execution"] = {"allow_unspliced_assembly": True}
    if consensus_mode:
        overrides["coding_potential"] = {"consensus": consensus_mode}
    try:
        cfg = load_config(cfg_file, overrides)
    except ConfigError as e:
        print(f"\nCannot run - please fix this and try again:\n   {e}")
        return

    if args.dry_run:
        print("\n(dry run) Configuration is valid - nothing was executed. "
              "Remove --dry-run to launch the run.")
        return

    # ---- STEP 1: make sure every required tool is installed ----
    if not setup_tools(cfg):
        print("\nStopping: not all required tools are available.")
        return

    # ---- STEP 2: make sure the genome + annotation are present ----
    if not setup_reference(cfg):
        print("\nStopping: genome / annotation not available (see above).")
        return

    # ---- STEP 3: run the pipeline (downloads the SRA data, then processes it) ----
    print("\n" + "-" * 52)
    print("STEP 3 - run the pipeline (reads -> alignment -> ... -> lncRNAs)")
    print("-" * 52)
    print("Everything is ready - starting the run.")
    # Save the resolved config so an interrupted run can be resumed without re-asking.
    save_run_config(cfg)
    _launch_engine(cfg)


if __name__ == "__main__":
    main()
