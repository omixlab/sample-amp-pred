# Benchmark results

`sample-amp-pred` scored on the external-benchmark cohort (n=967, 126 positive, 13.0% prevalence), at threshold 0.5:

| Metric | Value |
|---|---|
| roc_auc | 0.9321 |
| pr_auc | 0.6947 |
| accuracy | 0.9028 |
| precision | 0.6290 |
| recall | 0.6190 |
| specificity | 0.9453 |
| f1 | 0.6240 |
| mcc | 0.5682 |

Confusion matrix at threshold 0.5: TP=78, FP=46, FN=48, TN=795.

This matches sAMPle's own published numbers for this same cohort/model
(accuracy/precision/recall/F1/MCC exactly, to 4 decimals) with one small, expected exception:
ROC-AUC and PR-AUC here are 0.932/0.695 against sAMPle's own published 0.933/0.696 — a ~0.001
difference. That is consistent with this package's one documented limitation (see the main
README and `CHECKS.md`): a handful of these benchmark peptides are disulfide-bridged or cyclic
and get a small molecular-weight/charge correction in sAMPle's own training-time features, which
this standalone tool does not reproduce. The correction is too small to flip any classification
at threshold 0.5 (hence the identical confusion matrix), but it nudges a few raw probabilities
enough to move the threshold-free ranking metrics by about a thousandth.

The input FASTA holds 968 sequences; one (`SAMPLE_DECOY_00049796`) carries the non-standard
residue `B` and is reported as `prediction=error` rather than scored, leaving exactly 967 scored
sequences — the same cohort size sAMPle's own manuscript reports for its model-vs-Macrel/
AMPlify/AMPidentifier comparison. (amPEPpy's own published comparison used a 968-sequence
cohort, one more than the other three tools, for its own compatibility reasons unrelated to this
one dropped sequence.)

## Context: sAMPle's own published comparison on this cohort

Reference tools' numbers below are sAMPle's own published aggregate metrics on the same cohort (no per-sequence scores for them are available here, so this is context, not a head-to-head run by this tool):

| Tool | ROC-AUC | PR-AUC |
|---|---|---|
| **sample-amp-pred (this run)** | **0.932** | **0.695** |
| Macrel | 0.887 | 0.585 |
| AMPlify | 0.907 | 0.607 |
| AMPidentifier | 0.881 | 0.549 |
| amPEPpy | 0.836 | 0.483 |
