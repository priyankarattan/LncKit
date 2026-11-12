"""The pipeline engine - orchestrates the whole run.

Responsibilities:
  * fetch SRA accessions (when the run starts from accessions) and discover samples,
  * run each upstream step through its configured backend,
  * extract novel transcripts by gffcompare class code,
  * chain the identification filters
        length -> exon -> expression -> ORF -> coding-potential -> protein-homology
        -> housekeeping ncRNA
    each narrowing one candidate set, recording in/out counts at every step,
  * classify the survivors and write the final lncRNA GTF/FASTA/ID list,
  * checkpoint + resume, and record everything into run metadata for the report.

Filter order is deliberate, and is a performance decision as much as a biological
one: the cheap, decisive filters (class code, length, exon count, expression, ORF)
run first and typically remove most of the assembly, so the expensive models (CPAT,
DIAMOND, cmscan) only ever see the survivors.

Execution model is "orchestrated pure-Python": explicit, debuggable, dependency-free,
with checkpoint-based resume. See DESIGN.md for why this was chosen over
Snakemake/Nextflow for the default engine (a Nextflow export is a planned add-on).
"""
from __future__ import annotations

import time
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Set

from lnckit.core.config import Config, compatibility_warnings
from lnckit.core.context import RunContext
from lnckit.core.logging_utils import get_logger, setup_logging
from lnckit.core.metadata import StepRecord
from lnckit.filters.exon import filter_by_exons
from lnckit.filters.expression import filter_by_expression
from lnckit.filters.length import filter_by_length
from lnckit.filters.orf import OrfConfig, filter_by_orf
from lnckit.utils.fasta import subset_fasta
from lnckit.utils.gffcompare import find_tmap, novel_transcript_ids
from lnckit.utils.gtf import (has_empty_transcript_ids, has_invalid_strand,
                              sanitize_gtf, subset_gtf_by_ids, transcript_ids)
from lnckit.workflow.steps import PIPELINE, STEP_NAMES

log = get_logger("lnckit.engine")


class Engine:
    def __init__(self, config: Config, dry_run: bool = False):
        self.cfg = config
        self.ctx = RunContext.from_config(config)
        self.dry_run = dry_run or config.get("execution.dry_run", False)

    # -- helpers ----------------------------------------------------------
    def _record(self, name, tool="", status="completed", t0=0.0,
                n_in=-1, n_out=-1, outputs=None):
        self.ctx.metadata.add_step(StepRecord(
            name=name, tool=tool or "", status=status,
            runtime_seconds=round(time.time() - t0, 2) if t0 else 0.0,
            transcripts_in=n_in, transcripts_out=n_out,
            outputs=[str(o) for o in (outputs or [])],
        ))

    def _apply(self, name: str, tool: str, result, candidates: Set[str],
               t0: float, criterion: str = "") -> Set[str]:
        """Record a FilterResult and return the surviving candidates.

        `criterion` is a short human phrase for what this filter keeps (e.g.
        "length >= 200 nt"); it is printed alongside the before/after sequence counts
        so each step reads like "5000 sequences in -> 4187 passed length >= 200 nt".

        `result is None` means the stage was skipped because the reference data it
        needs was not configured. That is recorded as 'skipped', never as a pass:
        the report must not imply a filter ran when it did not.
        """
        if result is None:
            log.info("%-16s skipped - %d sequence(s) carried through unchanged",
                     name, len(candidates))
            self._record(name, tool, status="skipped", t0=t0,
                         n_in=len(candidates), n_out=len(candidates))
            self._write_step_output(name, candidates)
            return candidates
        self._record(name, tool, t0=t0, n_in=result.n_in, n_out=result.n_out)
        crit = f" [keep {criterion}]" if criterion else ""
        log.info("%-16s %d sequence(s) in -> %d passed, %d removed%s",
                 name, result.n_in, result.n_out, len(result.removed), crit)
        self._write_step_output(name, result.kept)
        return result.kept

    def _write_step_output(self, step: str, kept: Set[str]) -> None:
        """Save the transcripts that PASSED this step, so every filtered file is on
        disk and visibly feeds the next step. Files land in filtered_transcripts/."""
        d = self.ctx.dirs["filtered_transcripts"]
        (d / f"{step}.ids.txt").write_text("\n".join(sorted(kept)) + "\n")
        if getattr(self, "_merged_gtf", None):
            subset_gtf_by_ids(self._merged_gtf, kept, d / f"{step}.gtf")
        if getattr(self, "_tx_fasta", None):
            subset_fasta(self._tx_fasta, kept, d / f"{step}.fasta")

    def plan(self) -> List[str]:
        return STEP_NAMES

    @contextmanager
    def _stage(self, name: str):
        """Log a step's start / completion / failure to console AND the run log.

        On failure it records the step as 'failed' in the run metadata and flushes it to
        disk before re-raising, so <output>/logs/lnckit.log and run_metadata.json always
        show exactly which step stopped the run and why.
        """
        log.info("[START]  %s", name)
        t0 = time.time()
        try:
            yield
        except Exception as e:
            log.error("[FAILED] %s - %s: %s", name, type(e).__name__, e)
            self._record(name, status="failed", t0=t0)
            try:
                self.ctx.metadata.write(self.ctx.output_dir)
            except Exception:               # never mask the real error with a write error
                pass
            raise
        log.info("[DONE]   %s  (%.1fs)", name, time.time() - t0)

    # -- main entry -------------------------------------------------------
    def run(self, from_step: Optional[str] = None, until_step: Optional[str] = None):
        # Enable the plain-text run log now that the output dir is known, so every step
        # start/done/failure line is captured to <output>/logs/lnckit.log as well.
        setup_logging(log_dir=self.ctx.logs,
                      level=self.cfg.get("logging.level", "INFO"))
        log.info("=== LncKit run started - output dir: %s ===", self.ctx.output_dir)

        for w in compatibility_warnings(self.cfg):
            log.warning(w)

        if self.dry_run:
            log.info("DRY RUN - planned steps:")
            for s in PIPELINE:
                log.info("  %-14s %s", s.name, s.description)
            return self.ctx

        names = STEP_NAMES
        start_i = names.index(from_step) if from_step else 0
        end_i = names.index(until_step) if until_step else len(names) - 1
        active = set(names[start_i:end_i + 1])

        self._validate_reference_files()
        self._sanitize_annotation()

        from lnckit.modules import align, assemble, compare, fetch, index, qc, trim

        # ---- inputs: SRA accessions -> FASTQ ----------------------------
        reads_dir = None
        if "fetch" in active:
            with self._stage("fetch"):
                t0 = time.time()
                reads_dir = fetch.run_fetch(self.ctx)
                self._record("fetch", self.cfg.get("sra.tool", "fastq-dump"), t0=t0)
        # Discover the PER-RUN samples (one per accession); QC and trimming operate on
        # each replicate individually - cleaner than trimming one giant merged file.
        samples = self._discover_samples(reads_dir)

        # ---- upstream ---------------------------------------------------
        if "qc" in active and samples:
            with self._stage("qc"):
                t0 = time.time()
                qc.run_qc(self.ctx, [s["fq1"] for s in samples])
                self._record("qc", self.cfg.tools.get("qc"), t0=t0)

        trimmed = samples
        if "trim" in active and samples:
            with self._stage("trim"):
                t0 = time.time()
                trimmed = trim.run_trim(self.ctx, samples)
                self._record("trim", self.cfg.tools.get("trim"), t0=t0)

        # ---- optionally pool replicates by group (cat) AFTER trimming --------------
        # Each group's TRIMMED runs are concatenated into one sample per group, so
        # assembly sees maximum depth per condition. Runs on the trimmed reads.
        align_input = trimmed
        if "merge" in active and self.cfg.get("inputs.merge_replicates", False):
            groups = self.cfg.get("inputs.groups", {}) or {}
            if groups and trimmed:
                from lnckit.modules import merge
                with self._stage("merge"):
                    t0 = time.time()
                    merged = merge.run_merge(
                        self.ctx, groups, trimmed, self.cfg.get("inputs.paired", True))
                    self._record("merge", "cat", t0=t0)
                if merged:
                    align_input = merged
            else:
                log.warning("inputs.merge_replicates set but no groups/trimmed reads "
                            "found - aligning per-run samples instead")

        if "index" in active:
            with self._stage("index"):
                index_ref = index.run_index(self.ctx)
        else:
            index_ref = None

        if "align" in active:
            with self._stage("align"):
                bams = align.run_align(self.ctx, align_input, index_ref)
        else:
            bams = self._existing_bams()

        # Resolve 'auto' strandedness from the aligned reads before assembly. Safe by
        # construction: detect() never raises and only returns RF/FR on overwhelming
        # evidence, otherwise None (unstranded) - so the worst case is today's safe default.
        if self.cfg.get("inputs.strandness") == "auto":
            try:
                from lnckit.modules import strandness as _sd
                probe = bams[0] if bams else None
                call, why = (_sd.detect(probe) if probe
                             else (None, "no BAM available - using unstranded"))
            except Exception as e:                 # never let detection stop the run
                call, why = None, f"detection skipped ({type(e).__name__}) - using unstranded"
            self.cfg.raw.setdefault("inputs", {})["strandness"] = call
            log.info("strandness auto-detect: %s", why)

        if "assemble" in active:
            with self._stage("assemble"):
                merged_gtf, tx_fasta = assemble.run_assemble(self.ctx, bams)
        else:
            merged_gtf, tx_fasta = self._existing_assembly()

        compare_prefix = None
        if "compare" in active and merged_gtf:
            with self._stage("compare"):
                t0 = time.time()
                compare_prefix = compare.run_compare(self.ctx, merged_gtf)
                self._record("compare", self.cfg.tools.get("compare", "gffcompare"), t0=t0)
        elif merged_gtf:
            prefix = self.ctx.tmp / "gffcompare" / "gffcmp"
            compare_prefix = prefix if find_tmap(prefix) else None

        if merged_gtf and tx_fasta:
            with self._stage("identify"):
                self._identify(merged_gtf, tx_fasta, compare_prefix, bams, active)

        self.ctx.metadata.write(self.ctx.output_dir)
        if "report" in active:
            with self._stage("report"):
                from lnckit.reports.report import build_report
                build_report(self.ctx.output_dir)
        log.info("=== Run complete. Outputs under %s ===", self.ctx.output_dir)
        return self.ctx

    def _feelnc_filter(self, merged_gtf: Path, candidates: Set[str]):
        """Run FEELnc_filter.pl (length + monoexonic + mRNA-overlap removal) as the
        candidate filter. Returns (surviving_ids, True) on success, or (candidates, False)
        to fall back to the built-in filters (FEELnc env missing / no annotation / error)."""
        import shutil
        from lnckit.wrappers import get_wrapper
        f = self.cfg.get("filters", {}) or {}
        annotation = self.cfg.get("reference.annotation_gtf")
        if not annotation or not shutil.which("FEELnc_filter.pl"):
            log.warning("candidate_filter=feelnc but FEELnc_filter.pl or reference "
                        "annotation is unavailable - using the built-in length/exon "
                        "filters instead.")
            return candidates, False
        min_len = int(f.get("min_transcript_length", 200))
        # LncKit min_exon_number: 1 keeps monoexonic -> FEELnc --monoex 1 (keep); >=2
        # multi-exon only -> --monoex 0.
        monoex = 1 if int(f.get("min_exon_number", 1)) <= 1 else 0
        out_gtf = self.ctx.dirs["filtered_transcripts"] / "feelnc_filter.gtf"
        try:
            with self._stage("candidate_filter"):
                t0 = time.time()
                get_wrapper("candidate_filter", "feelnc")(
                    threads=self.cfg.threads, log_dir=self.ctx.logs).run(
                    candidate_gtf=merged_gtf, annotation=annotation,
                    min_size=min_len, monoex=monoex, stdout_path=out_gtf)
                kept = candidates & transcript_ids(out_gtf)
                self._record("candidate_filter", "FEELnc_filter", t0=t0,
                             n_in=len(candidates), n_out=len(kept))
                log.info("%-16s %d sequence(s) in -> %d passed  "
                         "[FEELnc_filter: length >= %d nt, monoex=%d, mRNA-overlap removed]",
                         "candidate_filter", len(candidates), len(kept), min_len, monoex)
                self._write_step_output("candidate_filter", kept)
        except Exception as e:                       # noqa: BLE001
            log.warning("FEELnc_filter failed (%s) - using the built-in filters instead", e)
            return candidates, False
        return kept, True

    # -- identification core ---------------------------------------------
    def _identify(self, merged_gtf: Path, tx_fasta: Path,
                  compare_prefix: Optional[Path], bams: List[Path], active: Set[str]):
        from lnckit.classification.classifier import classify
        from lnckit.modules import identify, quantify

        # remember these so _apply() can save each step's filtered file
        self._merged_gtf = merged_gtf
        self._tx_fasta = tx_fasta

        f = self.cfg.get("filters", {}) or {}
        candidates = transcript_ids(merged_gtf)
        log.info("assembled transcripts: %d", len(candidates))

        # --- candidate filtering: FEELnc_filter OR the built-in cascade ---------------
        # tools.candidate_filter=feelnc runs FEELnc_filter.pl (length -s + monoexonic
        # --monoex + removal of transcripts overlapping known mRNA) in ONE step - the
        # published-pipeline approach. It REPLACES the built-in novel/length/exon steps.
        # Default 'builtin' uses gffcompare-novelty + filter_length + filter_exon.
        used_feelnc_filter = False
        if (self.cfg.tools.get("candidate_filter") == "feelnc"
                and "filter_length" in active):
            candidates, used_feelnc_filter = self._feelnc_filter(merged_gtf, candidates)

        # --- novel candidates by gffcompare class code -------------------
        if not used_feelnc_filter and compare_prefix:
            t0 = time.time()
            tmap = find_tmap(Path(compare_prefix))
            if tmap:
                keep_classes = f.get("keep_classes", ["u", "x", "i", "o"])
                kept = candidates & novel_transcript_ids(tmap, keep_classes)
                self._record("novel", "gffcompare", t0=t0,
                             n_in=len(candidates), n_out=len(kept))
                log.info("%-16s %7d -> %7d   (class codes %s)",
                         "novel", len(candidates), len(kept), sorted(keep_classes))
                candidates = kept
            else:
                log.warning("no gffcompare .tmap found beside %s - every assembled "
                            "transcript is treated as a candidate", compare_prefix)

        # --- cheap structural filters (skipped if FEELnc_filter already did them) ----
        if "filter_length" in active and not used_feelnc_filter:
            t0 = time.time()
            min_len = f.get("min_transcript_length", 200)
            r = filter_by_length(tx_fasta, candidates, min_len)
            candidates = self._apply("filter_length", "builtin", r, candidates, t0,
                                     criterion=f"length >= {min_len} nt")

        if "filter_exon" in active and not used_feelnc_filter:
            t0 = time.time()
            min_exon = f.get("min_exon_number", 1)
            r = filter_by_exons(merged_gtf, candidates, min_exon)
            candidates = self._apply("filter_exon", "builtin", r, candidates, t0,
                                     criterion=f">= {min_exon} exon(s)")

        # Expression metric: FPKM if min_expression_fpkm set (or expression_metric=fpkm),
        # else TPM if min_expression_tpm set. If NEITHER threshold is set the filter is
        # DISABLED (the default) - abundance does not define a lncRNA, so we skip it.
        fpkm_cut = f.get("min_expression_fpkm")
        tpm_cut = f.get("min_expression_tpm")
        want_fpkm = fpkm_cut is not None or str(f.get("expression_metric", "")).lower() == "fpkm"
        expr_enabled = want_fpkm or tpm_cut is not None
        if "filter_expr" in active and bams and expr_enabled:
            t0 = time.time()
            if want_fpkm:
                metric, min_expr = "fpkm", (fpkm_cut if fpkm_cut is not None else 1.0)
            else:
                metric, min_expr = "tpm", tpm_cut
            min_samp = f.get("min_expressed_samples", 1)
            abund = quantify.run_quantify(self.ctx, bams, merged_gtf, metric=metric)
            r = filter_by_expression(abund, candidates,
                                     min_tpm=min_expr, min_samples=min_samp)
            log.info("expression filter: %s >= %s in >= %d sample(s)",
                     metric.upper(), min_expr, min_samp)
            candidates = self._apply("filter_expr", f"stringtie/{metric}", r, candidates, t0,
                                     criterion=f"{metric.upper()} >= {min_expr} in >= "
                                               f"{min_samp} sample(s)")
        elif "filter_expr" in active and not expr_enabled:
            log.info("expression filter: disabled (no min_expression_fpkm/tpm set) - skipping")

        # --- ORF filter (model-independent coding gate) -------------------
        # DISABLED by default: the coding-potential tools already weigh ORF length/coverage,
        # so this separate getorf gate is redundant unless explicitly enabled (orf.enabled).
        orf_conf = f.get("orf", {}) or {}
        if "filter_orf" in active and not orf_conf.get("enabled", False):
            log.info("ORF filter: disabled (filters.orf.enabled=false) - the coding-potential "
                     "tools already assess ORFs; skipping getorf")
        if "filter_orf" in active and orf_conf.get("enabled", False):
            t0 = time.time()
            orf_cfg = OrfConfig(**{k: v for k, v in orf_conf.items()
                                   if k in OrfConfig.__annotations__})
            getorf_out, orf_tool = self._maybe_run_orf_tool(tx_fasta, candidates, orf_conf)
            r = filter_by_orf(tx_fasta, candidates, orf_cfg, getorf_output=getorf_out)
            # `orf_tool` is what actually ran, which is not always what was
            # configured: the metadata must not credit getorf for a fallback run.
            candidates = self._apply("filter_orf", orf_tool, r, candidates, t0,
                                     criterion="no long/complete ORF (non-coding)")

        # --- evidence stages ---------------------------------------------
        if "coding" in active:
            t0 = time.time()
            cp = self.cfg.tools.get("coding_potential")
            cp_label = "+".join(cp) if isinstance(cp, (list, tuple)) else cp
            r = identify.run_coding_potential(self.ctx, tx_fasta, candidates)
            candidates = self._apply("coding", cp_label, r, candidates, t0,
                                     criterion=f"non-coding by {cp_label}")

        if "protein" in active:
            t0 = time.time()
            r = identify.run_protein_filter(self.ctx, tx_fasta, candidates)
            candidates = self._apply("protein", self.cfg.tools.get("protein_filter"),
                                     r, candidates, t0,
                                     criterion="no protein-database match")

        if "housekeeping" in active:
            t0 = time.time()
            r = identify.run_housekeeping(self.ctx, tx_fasta, candidates)
            candidates = self._apply("housekeeping", self.cfg.tools.get("housekeeping"),
                                     r, candidates, t0,
                                     criterion="not a known housekeeping ncRNA")

        # --- final set ----------------------------------------------------
        final_dir = self.ctx.dirs["final"]
        final_gtf = final_dir / "lncRNA.gtf"
        final_fa = final_dir / "lncRNA.fasta"
        subset_gtf_by_ids(merged_gtf, candidates, final_gtf)
        subset_fasta(tx_fasta, candidates, final_fa)
        (final_dir / "lncRNA_ids.txt").write_text("\n".join(sorted(candidates)) + "\n")

        self.ctx.metadata.final_counts["high_confidence_lncRNA"] = len(candidates)
        self._record("identify", "chain", n_out=len(candidates),
                     outputs=[final_gtf, final_fa])
        log.info("high-confidence lncRNAs: %d -> %s", len(candidates), final_gtf)

        if "classify" in active and self.cfg.get("reference.annotation_gtf"):
            t0 = time.time()
            out = classify(self.ctx, final_gtf, compare_prefix)
            self._record("classify", self.cfg.tools.get("classifier"), t0=t0,
                         n_out=len(candidates), outputs=[out])
        return final_gtf

    def _validate_reference_files(self) -> None:
        """Catch mixed-up reference paths early with a clear message.

        A frequent mistake is pointing reference.genome_fasta at the sample sheet, a GTF,
        or a protein file. gffread/HISAT2 then fail with a cryptic error deep in the run
        ('sequence lines ... must have the same length'). We check the genome is really a
        nucleotide FASTA (first non-blank line starts with '>') before anything uses it.
        """
        genome = self.cfg.get("reference.genome_fasta")
        if not genome:
            return
        p = Path(genome)
        if not p.exists():
            return                                   # 'missing' is reported elsewhere
        try:
            with open(p) as fh:
                for line in fh:
                    if line.strip():
                        first = line.rstrip("\n")
                        break
                else:
                    first = ""
        except OSError:
            return
        if not first.startswith(">"):
            raise ConfigError(
                f"reference.genome_fasta = {genome}\n  This is not a FASTA - its first "
                f"line is:\n    {first[:80]!r}\n  You likely pointed it at your sample "
                f"sheet, a GTF/GFF, or a protein file by mistake. Set genome_fasta to the "
                f"genome NUCLEOTIDE FASTA (e.g. YourGenome.fasta, first line '>chr...').")

    @staticmethod
    def _looks_like_gff(ann: Path) -> bool:
        """True if `ann` is GFF/GFF3 rather than GTF (StringTie wants proper GTF).

        Decides by extension first, then by sniffing the attribute column: GTF uses
        `gene_id "x"; transcript_id "y";`, GFF3 uses `ID=x;Parent=y`.
        """
        suf = ann.suffix.lower()
        if suf in (".gff", ".gff3"):
            return True
        if suf == ".gtf":
            return False
        try:
            with open(ann) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9:
                        continue
                    attr = cols[8]
                    if 'transcript_id "' in attr or 'gene_id "' in attr:
                        return False            # GTF
                    return "=" in attr          # GFF3-style attributes
        except OSError:
            pass
        return False

    def _ensure_gtf(self, ann: Path) -> Path:
        """Return a GTF for `ann`, converting a GFF/GFF3 with gffread if needed.

        StringTie frequently can't extract reference transcripts from a GFF3
        ('could not find any valid reference transcripts'); gffread -T rewrites it into a
        clean GTF (proper gene_id/transcript_id + exon features) that StringTie reads.
        """
        if not self._looks_like_gff(ann):
            return ann
        from lnckit.utils.shell import run as run_tool
        gtf = self.ctx.output_dir / ".lnckit" / "annotation.from_gff.gtf"
        gtf.parent.mkdir(parents=True, exist_ok=True)
        resume = self.cfg.get("execution.resume", True)
        fresh = (gtf.exists() and resume and gtf.stat().st_size > 0
                 and gtf.stat().st_mtime >= ann.stat().st_mtime)
        if not fresh:
            log.warning("annotation %s is GFF/GFF3; converting to GTF with gffread "
                        "(StringTie needs GTF) -> %s", ann, gtf)
            run_tool(["gffread", str(ann), "-T", "-o", str(gtf)], tool="gffread",
                     log_file=self.ctx.logs / "gffread.annotation.log")
        if not gtf.exists() or gtf.stat().st_size == 0:
            raise ConfigError(
                f"Converting {ann} to GTF produced nothing - the file has no usable "
                f"transcript/exon features. Check it is a real coordinate GFF3 (columns: "
                f"seqid, source, type[gene/mRNA/exon/CDS], start, end, ...), not a FASTA "
                f"or a features-only file. Inspect with: cut -f3 '{ann}' | sort | uniq -c")
        return gtf

    def _sanitize_annotation(self) -> None:
        """Repoint reference.annotation_gtf at a StringTie-safe copy when needed.

        Two problems are handled: (1) a GFF/GFF3 annotation is converted to GTF with
        gffread (StringTie can't reliably read GFF3); (2) NCBI RefSeq GTFs carry `gene`
        rows with an empty `transcript_id ""`, which makes StringTie abort - those are
        cleaned. The result is written under .lnckit/ and used for the rest of the run
        (assembly, merge, gffcompare, classification). No-op for an already-clean GTF.
        """
        ann = self.cfg.get("reference.annotation_gtf") or \
            self.cfg.get("reference.annotation_gff")
        if not ann:
            return
        ann = Path(ann)
        if not ann.exists():
            return
        ann = self._ensure_gtf(ann)                       # GFF/GFF3 -> GTF if needed
        if has_empty_transcript_ids(ann) or has_invalid_strand(ann):
            clean = self.ctx.output_dir / ".lnckit" / "annotation.clean.gtf"
            clean.parent.mkdir(parents=True, exist_ok=True)
            resume = self.cfg.get("execution.resume", True)
            # Reuse a cached clean copy only if it is newer than the source AND itself
            # already clean - a copy left by an older sanitizer may still carry the
            # rows we now strip (e.g. '?'-strand rps12), so regenerate it in that case.
            fresh = (clean.exists() and resume
                     and clean.stat().st_mtime >= ann.stat().st_mtime
                     and not has_invalid_strand(clean)
                     and not has_empty_transcript_ids(clean))
            if not fresh:
                c = sanitize_gtf(ann, clean)
                log.warning("annotation %s: dropped %d row(s) StringTie's -G reader "
                            "rejects (%d empty transcript_id, %d undetermined-strand '?', "
                            "e.g. trans-spliced chloroplast rps12); using cleaned copy %s",
                            ann, c["total"], c["no_tid"], c["bad_strand"], clean)
            else:
                log.info("reusing cleaned annotation from a previous run: %s", clean)
            ann = clean
        self.cfg.raw.setdefault("reference", {})["annotation_gtf"] = str(ann)

    def _maybe_run_orf_tool(self, tx_fasta: Path, candidates: Set[str],
                            orf_conf: dict):
        """Run the external ORF caller if available; else fall back to pure Python.

        Returns (getorf_output_or_None, tool_that_actually_ran) so the run metadata
        records the ORF source truthfully rather than the one that was requested.
        """
        if not orf_conf or not orf_conf.get("enabled", True):
            return None, "disabled"
        tool = orf_conf.get("tool", "getorf")
        if tool not in {"getorf", "orffinder"}:
            # TransDecoder emits peptides, not the coordinate headers the parser
            # reads, so the built-in finder is used instead of mis-parsing them.
            return None, "python-builtin"
        try:
            from lnckit.wrappers import get_wrapper
            fasta = subset_fasta(tx_fasta, candidates,
                                 self.ctx.dirs["orf"] / "orf_input.fasta")
            out = self.ctx.dirs["orf"] / "orfs.fa"
            get_wrapper("orf", tool)(threads=self.cfg.threads, log_dir=self.ctx.logs).run(
                fasta=fasta, orf_out=out, min_nt=30,
                genetic_code=orf_conf.get("genetic_code", 1),
                strand=orf_conf.get("strand", "both"))
            return out, tool
        except Exception as e:      # tool missing/failed -> transparent fallback
            log.warning("ORF caller '%s' unavailable (%s); using the built-in Python "
                        "ORF finder instead.", tool, e)
        return None, "python-builtin"

    # -- artefacts from a previous run (resume / single-stage commands) -----
    def _existing_bams(self) -> List[Path]:
        """BAMs already on disk, so `lnckit assemble` or `resume` can pick up."""
        bams = sorted(self.ctx.dirs["alignment"].glob("*.bam"))
        if bams:
            log.info("reusing %d BAM(s) from a previous run", len(bams))
        return bams

    def _existing_assembly(self):
        """(merged_gtf, transcript_fasta) from a previous run, if both are present."""
        merged = self.ctx.dirs["assembly"] / "merged.gtf"
        fasta = self.ctx.dirs["assembly"] / "transcripts.fasta"
        if merged.exists() and fasta.exists():
            log.info("reusing assembly from a previous run: %s", merged)
            return merged, fasta
        return None, None

    # -- input discovery --------------------------------------------------
    def _discover_samples(self, reads_dir: Optional[Path] = None) -> List[Dict]:
        rd = reads_dir or self.cfg.get("inputs.reads_dir")
        if not rd:
            log.warning("Neither inputs.reads_dir nor inputs.sra_accessions is set; "
                        "there is nothing to process.")
            return []
        rd = Path(rd)
        paired = self.cfg.get("inputs.paired", True)
        samples: List[Dict] = []
        if paired:
            for r1 in sorted(rd.glob("*_1.f*q*")):
                name = r1.name.split("_1.")[0]
                r2 = next(rd.glob(f"{name}_2.f*q*"), None)
                if r2 is None:
                    log.warning("%s has no mate file - skipped (set inputs.paired: false "
                                "for single-end data)", r1.name)
                    continue
                samples.append({"name": name, "fq1": r1, "fq2": r2})
        else:
            for r1 in sorted(rd.glob("*.f*q*")):
                if "_2." in r1.name:
                    continue
                samples.append({"name": r1.name.split(".f")[0], "fq1": r1, "fq2": None})
        log.info("discovered %d sample(s) in %s", len(samples), rd)
        return samples
