"""LncFinder wrapper - SVM coding-potential classifier (R / Bioconductor).

LncFinder is an R package, not a command-line tool, so there is no binary to call.
The wrapper writes a small R driver script and runs it with `Rscript`; the science
still lives entirely in the LncFinder package - lnckit only marshals input/output.

Plant use: LncFinder's built-in models were trained on animal data, so for plants
supply the retrained LncFinder-plant model (and its frequencies file) via
`coding_potential.lncfinder.model` / `.frequencies` - mirroring how CPAT takes the
plant hexamer/logit model. Without them the wrapper falls back to a built-in model
and warns that it is not plant-tuned.
"""
from __future__ import annotations
import shutil
from functools import lru_cache
from pathlib import Path
from typing import List, Optional

from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@lru_cache(maxsize=None)
def _lncfinder_rscript() -> str:
    """Rscript whose R library actually has LncFinder.

    LncFinder + its heavy R deps live in a dedicated `lnckit_lncfinder` conda env
    (kept apart so a solve never perturbs the main tool env). It sits beside the other
    envs, so we locate it from any Rscript on PATH and prefer it; if that env is absent
    we fall back to the plain `Rscript` (the run then simply skips LncFinder).
    """
    main = shutil.which("Rscript")
    if main:
        cand = Path(main).parents[2] / "lnckit_lncfinder" / "bin" / "Rscript"
        if cand.exists():
            return str(cand)
    return main or "Rscript"

# R driver: read FASTA -> lnc_finder() -> write ID/Pred/Coding_prob TSV.
# model : a built-in name ("mouse") OR a path to a saved model. Plant-LncPipe's
#         Plant_model.rda is saved with saveRDS, so we try readRDS() FIRST, then fall
#         back to load() for classic .RData objects - either works transparently.
# freq  : a saved frequencies object (path) OR the special form "make:mRNA.fa,lncRNA.fa"
#         which recomputes frequencies with make_frequencies() from training FASTAs,
#         exactly like Plant-LncPipe does. Defaults to `model` if not given.
_R_DRIVER = r'''#!/usr/bin/env Rscript
suppressMessages(library(LncFinder))
suppressMessages(library(seqinr))
args  <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]; out <- args[2]
model <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "mouse"
freq  <- if (length(args) >= 4 && nzchar(args[4])) args[4] else model
cores <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1L

read_model <- function(x) {
  if (!file.exists(x)) return(x)                       # a built-in name, e.g. "mouse"
  obj <- tryCatch(readRDS(x), error = function(e) NULL) # saveRDS form (Plant_model.rda)
  if (!is.null(obj)) return(obj)
  e <- new.env(); nm <- load(x, envir = e); get(nm[1], envir = e)   # classic .RData
}
get_freq <- function(x, model_obj) {
  if (startsWith(x, "make:")) {                        # recompute from training FASTAs
    parts <- strsplit(sub("^make:", "", x), ",", fixed = TRUE)[[1]]
    mrna <- seqinr::read.fasta(trimws(parts[1]))
    lnc  <- seqinr::read.fasta(trimws(parts[2]))
    return(make_frequencies(cds.seq = mrna, lncRNA.seq = lnc, SS.features = FALSE,
                            cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE,
                            ignore.illegal = TRUE))
  }
  if (file.exists(x)) return(read_model(x))            # a saved frequencies object
  model_obj                                            # last resort: reuse the model name
}
mdl  <- read_model(model)
frq  <- get_freq(freq, mdl)
seqs <- seqinr::read.fasta(fasta)
res  <- lnc_finder(seqs, SS.features = FALSE, format = "DNA",
                   frequencies.file = frq, svm.model = mdl, parallel.cores = cores)
df <- data.frame(ID = rownames(res), Pred = res[["Pred"]],
                 Coding_prob = res[["Coding.Potential"]],
                 stringsAsFactors = FALSE)
write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
'''


@register("coding_potential", "lncfinder")
class LncFinderWrapper(BaseWrapper):
    requires = ("lncfinder",)          # resolves to Rscript; see dependencies.KNOWN_TOOLS
    tool_name = "lncfinder"

    def build_command(self, *, fasta: Path, out_tsv: Path,
                      model: str = "mouse", frequencies: Optional[str] = None,
                      **_) -> List[str]:
        script = Path(out_tsv).parent / "lncfinder_predict.R"
        script.parent.mkdir(parents=True, exist_ok=True)
        script.write_text(_R_DRIVER)
        return [_lncfinder_rscript(), str(script), str(fasta), str(out_tsv),
                str(model or ""), str(frequencies or model or ""), self.threads]

    def expected_outputs(self, *, out_tsv: Path, **_) -> List[Path]:
        return [Path(out_tsv)]
