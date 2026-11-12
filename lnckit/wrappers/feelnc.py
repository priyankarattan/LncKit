"""FEELnc wrapper - candidate filtering and genomic-location classification.

FEELnc_filter removes short/mono-exonic transcripts and those overlapping known
mRNAs; FEELnc_classifier assigns lncRNA classes relative to reference genes.

FEELnc lives in its own conda env (``lnckit_feelnc``): it pins an old Perl 5.22 with
its own bundled ``Bio::SeqFeature`` modules, and lnckit.sh puts that env on PATH *after*
the main ``lnckit`` env so its old Perl never shadows anything else. The catch is that
a bare ``FEELnc_*.pl`` (shebang ``#!/usr/bin/env perl``) would then run under the main
env's newer Perl, which lacks FEELnc's modules. So every FEELnc script is launched
explicitly with its sibling ``perl`` and a ``PERL5LIB`` pointing at that env's
site_perl - immune to PATH order and to a polluted inherited ``PERL5LIB``.
"""
from __future__ import annotations
import shutil
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple
from lnckit.wrappers.base import BaseWrapper
from lnckit.wrappers.registry import register


@lru_cache(maxsize=None)
def _feelnc_runtime(script: str) -> Tuple[Tuple[str, ...], Optional[str]]:
    """Return (argv_prefix, PERL5LIB) to run a FEELnc script under its own Perl.

    Resolves ``script`` on PATH, then runs it explicitly via the sibling ``perl`` with
    PERL5LIB set to the env's site_perl (located by finding FEELnc's own bundled module).
    Falls back to ``((script,), None)`` - the plain invocation - when the env cannot be
    located, so nothing breaks outside the lnckit.sh layout.
    """
    path = shutil.which(script)
    if not path:
        return (script,), None
    env_bin = Path(path).parent
    env_root = env_bin.parent
    perl = env_bin / "perl"
    mods = sorted(env_root.glob("lib/perl5/site_perl/*/Bio/SeqFeature/database_part.pm"))
    perl5lib = str(mods[0].parents[2]) if mods else None   # .../site_perl/<ver>
    prefix = (str(perl), path) if perl.exists() else (path,)
    return prefix, perl5lib


@register("candidate_filter", "feelnc")
class FeelncFilterWrapper(BaseWrapper):
    requires = ("feelnc",)
    tool_name = "FEELnc_filter"

    def build_command(self, *, candidate_gtf: Path, annotation: Path,
                      min_size: int = 200, monoex: int = -1, **_) -> List[str]:
        # writes result to stdout -> the step redirects it to a file.
        prefix, _ = _feelnc_runtime("FEELnc_filter.pl")
        return [*prefix, "-i", str(candidate_gtf), "-a", str(annotation),
                "--monoex", monoex, "-s", min_size, "-p", self.threads]

    def env(self, **_) -> Optional[dict]:
        _, perl5lib = _feelnc_runtime("FEELnc_filter.pl")
        return {"PERL5LIB": perl5lib} if perl5lib else None

    def expected_outputs(self, **_) -> List[Path]:
        return []  # stdout captured by caller


@register("classifier", "feelnc")
class FeelncClassifierWrapper(BaseWrapper):
    requires = ("feelnc",)
    tool_name = "FEELnc_classifier"

    def build_command(self, *, lncrna_gtf: Path, annotation: Path, **_) -> List[str]:
        prefix, _ = _feelnc_runtime("FEELnc_classifier.pl")
        return [*prefix, "-i", str(lncrna_gtf), "-a", str(annotation)]

    def env(self, **_) -> Optional[dict]:
        _, perl5lib = _feelnc_runtime("FEELnc_classifier.pl")
        return {"PERL5LIB": perl5lib} if perl5lib else None

    def expected_outputs(self, **_) -> List[Path]:
        return []
