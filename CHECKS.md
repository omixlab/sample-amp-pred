# Verification record

This tool's feature engineering (`src/sample_amp_pred/features.py`) is a from-scratch,
dependency-light reproduction of sAMPle's own `sample/inference.py`. It was checked against
sAMPle's own live inference path, not trusted by inspection alone.

## What was checked

Ran both implementations — this package's `features.build_feature_matrix` +
`bundle.estimator.predict_proba`, and sAMPle's own `sample.inference.build_feature_frame` +
`sample.model_registry.load_model`'s estimator — on the same set of sequences, and compared the
resulting probabilities.

Sequences checked: six real peptides, several drawn from the model's own holdout FASTA
(`flaml.holdout.fasta`), plus magainin-2 (a well-known, independently published AMP not
necessarily in that holdout) and a 2-residue minimum-length edge case.

```
sequence                                        real_proba   mine_proba   abs_diff
GNNRPVYIPQPRPPHPRI                              0.96524560   0.96524560   0.00e+00
GLFDIVKKVVGTLAGL                                0.97248060   0.97248060   0.00e+00
KSSAYSLQMGATAIKQVKKLFKKWGW                      0.85016578   0.85016578   0.00e+00
GIGASILSAGKSALKGLAKGLAEHFAN                     0.86659837   0.86659837   0.00e+00
GIGTKILGGVKTALKGALKELASTYAN                     0.89744908   0.89744908   0.00e+00
KWKLFKKIGIGAVLKVLTTGLPALIS                      0.98232263   0.98232263   0.00e+00
GIGKFLHSAKKFGKAFVGEIMNS                         0.90583098   0.90583098   0.00e+00
AA                                              0.71518952   0.71518952   0.00e+00

max abs diff across 8 sequences: 0.00e+00
```

Every probability matched to the full float precision printed (bitwise-identical, not merely
close). This also confirms:

- the bundled model's column order (`bundle.columns`) is identical to the live model's own
  (`sample.model_registry.load_model`'s second return value) — checked with a direct tuple
  equality assertion before scoring anything, so a silent column-order drift would have failed
  loudly rather than producing a quietly wrong number;
- the exported artifact's raw `xgboost.XGBClassifier` (`bundle["model"]._estimator` in sAMPle's
  own pickle) predicts identically to scoring through sAMPle's full `FlamlXGBoostEstimator`
  wrapper, confirming the wrapper is a pass-through with no extra transform this tool needs to
  reproduce.

## What was not checked (and why that's fine for this tool's scope)

- **Cyclic / disulfide-bridged peptides in sAMPle's own curated topology index.** sAMPle
  corrects molecular weight and net charge for those specific, previously-curated sequences.
  That correction depends on an internal database not published alongside the model, so it is
  not reproduced here — every sequence is scored as a linear chain, which is also what sAMPle's
  own server does for any sequence *outside* that index (i.e., every genuinely new sequence a
  user of this tool is likely to submit).
- **The applicability-domain / in-domain flag.** sAMPle's server reports whether a query sits
  inside the training data's feature-space envelope. This tool does not compute or report it,
  since the user-facing contract (`-i`/`-o`/`-f`/`-t`) asks only for probability and label.

## How to re-run this check

The cross-check script needs the full sAMPle repository on `sys.path` (to import
`sample.inference` and `sample.model_registry`) and is therefore not shipped in this
package — it is a development-time check, run once against sAMPle's own
`data/sample_final_release/models/antibacterial/general/any/classification/` artifacts at
export time, not something this package's own CI needs to repeat.
