"""Layered configuration system (YAML/JSON) with validation + compatibility checks.

Resolution order (lowest -> highest precedence):
    packaged defaults  <  user config file  <  CLI overrides
Everything is merged into one Config object that is passed, read-only, through
the whole run. The same object is serialised into run metadata for reproducibility.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from lnckit.core.exceptions import CompatibilityError, ConfigError

_PKG_DEFAULT = Path(__file__).parent.parent / "config" / "default.yaml"
_COMPAT = Path(__file__).parent.parent / "config" / "compatibility.yaml"

# Allowed choices per stage - the CLI/validation layer enforces these.
CHOICES = {
    "qc":             ["fastqc", "fastp", "multiqc"],
    "trim":           ["fastp", "trimmomatic", "cutadapt"],
    "aligner":        ["hisat2", "star", "bowtie2", "bwa"],
    "assembler":      ["stringtie", "cufflinks"],
    "coding_potential": ["cpat", "cpc2", "cnci", "plek", "lncfinder", "rnaplonc"],
    "protein_filter": ["diamond", "blastx"],
    "housekeeping":   ["infernal", "blastn", "none"],
    "classifier":     ["feelnc", "gffcompare"],
}

# Splice-aware aligners are required for reference-guided transcript assembly.
SPLICE_AWARE = {"hisat2", "star"}


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


@dataclass
class Config:
    """Validated, merged run configuration."""
    raw: Dict[str, Any] = field(default_factory=dict)

    # ---- convenience accessors -------------------------------------------
    def get(self, dotted: str, default=None):
        node = self.raw
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    @property
    def tools(self) -> Dict[str, str]:
        return self.raw.get("tools", {})

    @property
    def threads(self) -> int:
        return int(self.raw.get("resources", {}).get("threads", 4))

    @property
    def output_dir(self) -> Path:
        return Path(self.raw.get("output_dir", "lnckit_output"))

    def params_hash(self) -> str:
        blob = json.dumps(self.raw, sort_keys=True).encode()
        return hashlib.sha1(blob).hexdigest()[:12]

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.raw, sort_keys=False)


def load_config(
    user_path: Optional[Path] = None,
    overrides: Optional[Dict[str, Any]] = None,
) -> Config:
    """Build the effective Config from defaults + file + CLI overrides, then validate."""
    with open(_PKG_DEFAULT) as fh:
        merged = yaml.safe_load(fh) or {}

    if user_path:
        p = Path(user_path)
        if not p.exists():
            raise ConfigError(f"Config file not found: {p}")
        text = p.read_text()
        try:
            user = json.loads(text) if p.suffix.lower() == ".json" else yaml.safe_load(text)
        except Exception as e:
            raise ConfigError(f"Could not parse config {p}: {e}") from e
        merged = _deep_merge(merged, user or {})

    if overrides:
        merged = _deep_merge(merged, overrides)

    # Highest-precedence organism override (from lnckit.sh --organism / LNCKIT_ORGANISM). Applied
    # on EVERY load - fresh, batch (no prompt), AND resume (saved config) - so `organism: plant`
    # reliably reaches the run and the plant CPAT/LncFinder/CNCI models are used.
    env_org = os.environ.get("LNCKIT_ORGANISM")
    if env_org:
        merged["organism"] = env_org.strip().lower()

    # Same for strandness (lnckit.sh --strandness / LNCKIT_STRANDNESS). Reliably sets
    # inputs.strandness on EVERY load so StringTie gets --rf/--fr even in batch/resume.
    # 'none'/'unstranded' -> null (StringTie infers strand from splice motifs only).
    env_strand = os.environ.get("LNCKIT_STRANDNESS")
    if env_strand:
        v = env_strand.strip().lower()
        merged.setdefault("inputs", {})["strandness"] = (
            None if v in ("none", "null", "unstranded", "") else
            "auto" if v == "auto" else env_strand.strip().upper())

    _normalize_tool_lists(merged)
    _expand_paths(merged)
    cfg = Config(raw=merged)
    validate(cfg)
    return cfg


def _normalize_tool_lists(merged: dict) -> None:
    """Forgive a common config mistake: the TOOL LIST goes under `tools.coding_potential`,
    but `coding_potential:` (top level) is the SETTINGS block (consensus, cutoffs). If a
    user writes `coding_potential: [cpat, cpc2, ...]` at the top level (natural guess),
    that would silently NOT change which tools run - the default full ensemble (incl. PLEK)
    would still execute. Detect that and route the list to tools.coding_potential.
    """
    cp = merged.get("coding_potential")
    if isinstance(cp, (list, tuple)):
        merged.setdefault("tools", {})["coding_potential"] = list(cp)
        merged["coding_potential"] = {}       # leave the settings block empty (defaults)
        import logging
        logging.getLogger("lnckit.config").warning(
            "coding_potential was given as a list at the TOP level; using it as the tool "
            "selection (tools.coding_potential = %s). Put it under `tools:` to avoid this "
            "notice.", list(cp))

    # Legacy migration: 'pfam' used to be an hmmscan coding-potential backend. It is now a
    # diamond BLASTX homology database instead, so silently drop it from coding_potential
    # rather than hard-failing configs that still list it.
    tools = merged.get("tools")
    if isinstance(tools, dict):
        cpl = tools.get("coding_potential")
        if isinstance(cpl, (list, tuple)) and "pfam" in cpl:
            tools["coding_potential"] = [c for c in cpl if c != "pfam"]


def _expand_paths(merged: dict) -> None:
    """Expand `~` and $ENV_VARS in path fields, so configs stay machine-independent.

    Nobody's absolute layout is baked in: a user can write `~/dbs/uniprot.fasta` or
    `$LNCKIT_DB_DIR/Rfam.cm` and it resolves to THEIR paths at load time. Only real
    path keys are touched (never the `*_url` download links).
    """
    def expand(v):
        if not isinstance(v, str):
            return v
        out = os.path.expanduser(os.path.expandvars(v))
        # An unset $VAR is left literal by expandvars; don't hand a broken path
        # downstream - treat it as unconfigured (None) so the step degrades cleanly.
        return None if "$" in out and "${" in v else out

    for k in ("output_dir", "tmp_dir"):
        if isinstance(merged.get(k), str):
            merged[k] = expand(merged[k])
    ref = merged.get("reference")
    if isinstance(ref, dict):
        for k, v in ref.items():
            if k.endswith("_url"):
                continue
            if isinstance(v, str):
                ref[k] = expand(v)
            elif isinstance(v, (list, tuple)):    # e.g. protein_db: [uniprot, cog, ...]
                ref[k] = [expand(x) if isinstance(x, str) else x for x in v]
    inp = merged.get("inputs")
    if isinstance(inp, dict):
        for k in ("reads_dir", "samples", "sra_accession_file"):
            if isinstance(inp.get(k), str):
                inp[k] = expand(inp[k])
    # LncFinder model/frequencies can be $LNCKIT_DB_DIR-relative too (the plant model +
    # training FASTAs live in the package databases/). The frequencies "make:a.fa,b.fa"
    # form has its ${VARS} expanded in place by expandvars.
    lf = (merged.get("coding_potential") or {}).get("lncfinder") \
        if isinstance(merged.get("coding_potential"), dict) else None
    if isinstance(lf, dict):
        for k in ("model", "frequencies"):
            if isinstance(lf.get(k), str):
                lf[k] = expand(lf[k])


#: stages that accept a LIST of tools (run all, intersect the results) as well as a
#: single tool. `coding_potential` does so to reproduce the Plant-LncPipe ensemble
#: (CPAT and LncFinder): a transcript survives only if every model calls it non-coding.
ENSEMBLE_STAGES = {"coding_potential"}


def validate(cfg: Config) -> None:
    """Validate choices and cross-tool compatibility. Raises ConfigError/CompatibilityError."""
    tools = cfg.tools
    for stage, allowed in CHOICES.items():
        chosen = tools.get(stage)
        if chosen is None:
            continue
        is_list = isinstance(chosen, (list, tuple))
        if is_list and stage not in ENSEMBLE_STAGES:
            raise ConfigError(f"'{stage}' takes a single tool, not a list: {chosen!r}")
        for c in (chosen if is_list else [chosen]):
            if c not in allowed:
                raise ConfigError(
                    f"Invalid choice for '{stage}': {c!r}. Allowed: {allowed}"
                )

    # Compatibility rule: reference-guided assembly needs a splice-aware aligner.
    aligner = tools.get("aligner")
    assembler = tools.get("assembler")
    if assembler in {"stringtie", "cufflinks"} and aligner not in SPLICE_AWARE:
        # The user can override this from the wizard ("continue anyway") when they
        # accept that a non-splice-aware aligner produces unreliable assembly.
        if not cfg.get("execution.allow_unspliced_assembly", False):
            raise CompatibilityError(
                f"Transcript assembly with '{assembler}' requires a splice-aware aligner "
                f"(one of {sorted(SPLICE_AWARE)}), but '{aligner}' was selected. "
                f"Switch the aligner to HISAT2 or STAR, or choose a different assembler."
            )

    _check_backends_implemented(tools)


# Stages whose every allowed choice must have a registered wrapper. Excluded:
# housekeeping ("none" is valid but unregistered) and classifier ("gffcompare" is
# handled inline, not via the registry) - those are special-cased in their modules.
_MUST_BE_REGISTERED = {
    "qc": "qc", "trim": "trim", "aligner": "aligner", "assembler": "assembler",
    "coding_potential": "coding_potential", "protein_filter": "protein_filter",
}

# Backends handled inline in their module rather than by a single-command wrapper, so
# they legitimately have no entry under their registry stage (skip the registered-check).
_INLINE_BACKENDS = {"pfam", "rnaplonc"}   # coding_potential backends orchestrated in-module
#   pfam     : getorf translate + hmmscan   |   rnaplonc: txCdsPredict + feature_extraction + weka


def _check_backends_implemented(tools: Dict[str, str]) -> None:
    """Fail at config time if a selected backend has no wrapper, not mid-run.

    A choice listed in CHOICES but never implemented would otherwise pass validation
    and then raise KeyError from get_wrapper() hours into the pipeline.
    """
    try:
        from lnckit.wrappers import available
    except Exception:  # pragma: no cover - never let a registry import break config
        return
    for stage, reg_stage in _MUST_BE_REGISTERED.items():
        chosen = tools.get(stage)
        if not chosen:
            continue
        for c in (chosen if isinstance(chosen, (list, tuple)) else [chosen]):
            if c in _INLINE_BACKENDS:
                continue
            if c not in available(reg_stage):
                raise ConfigError(
                    f"'{stage}: {c}' is selected but no backend is registered for it. "
                    f"Available: {sorted(available(reg_stage))}."
                )


def compatibility_warnings(cfg: Config) -> List[str]:
    """Non-fatal advisories (e.g. suboptimal but allowed combinations)."""
    warns: List[str] = []
    tools = cfg.tools
    if tools.get("aligner") == "star" and cfg.get("resources.memory_gb", 32) < 30:
        warns.append("STAR typically needs >=30 GB RAM for plant genomes; "
                     "consider HISAT2 for low-memory nodes.")
    if tools.get("housekeeping") == "none":
        warns.append("Housekeeping ncRNA removal is disabled; rRNA/tRNA/snoRNA "
                     "contamination may remain in the final set.")
    return warns
