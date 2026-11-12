# Vendored third-party tools

## CPC2_standalone/
CPC2 (Coding Potential Calculator 2), from https://github.com/gao-lab/CPC2_standalone,
vendored so LncKit works offline (no clone at run time). It is **Python 2** and runs
under the `lnckit_py2` conda env (biopython + numpy). `setup.sh` compiles its bundled
libsvm (`libs/libsvm/libsvm-3.18.tar.gz`) and installs a `CPC2.py` shim that runs it
under that env. CPC2 is **GPL-3** - see `CPC2_standalone/LICENSE`. It is an independent
aggregated component; LncKit itself remains MIT.

## RNAplonc/
RNAplonc (plant lncRNA REPTree/Weka classifier), from
https://github.com/TatianneNegri/RNAplonc, vendored so it works offline. Chain:
txCdsPredict (bundled x86-64 binary) -> feature_extraction.pl (needs BioPerl, reused from
the lnckit_feelnc env) -> weka.jar REPTree with RNAplonc.model. Needs java (openjdk, in the
lnckit env). GPL-3 - see CPC2/RNAplonc upstream licenses. Independent aggregated component;
LncKit itself remains MIT.
