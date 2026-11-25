# Benchmark artifacts

This directory holds the external-benchmark cohort sAMPle's own manuscript scores its general
antibacterial classifier against (alongside Macrel, AMPlify, AMPidentifier and amPEPpy), the
predictions `sample-amp-pred` produces on it, and the resulting metrics and plots — so the
numbers in the main `README.md` can be checked against something concrete rather than taken on
faith.

| File | What it is |
|---|---|
| `external_holdout.fasta` | The benchmark cohort: 968 peptide sequences, a 40%-CD-HIT-partitioned sample of the model's own training holdout. |
| `external_holdout_labels.csv` | Ground-truth label (1 = active, 0 = inactive/decoy) for each sequence above. |
| `predictions.csv` | `sample-amp-pred`'s own output on `external_holdout.fasta` (`sample-amp-pred -i external_holdout.fasta -o predictions.csv`). |
| `reference_tool_metrics.json` | sAMPle's own published aggregate metrics for Macrel/AMPlify/AMPidentifier/amPEPpy on this same cohort — context only; no per-sequence scores for those tools are available here, so this package cannot re-run that part. |
| `results.json` / `results.md` | Metrics computed from `predictions.csv` against `external_holdout_labels.csv` — this package's own fresh numbers, not copied from the manuscript. |
| `plots/roc_curve.png`, `plots/pr_curve.png` | ROC and Precision-Recall curves from `predictions.csv`. |
| `plots/tool_comparison.png` | `sample-amp-pred`'s measured ROC-AUC/PR-AUC alongside the four reference tools' published numbers. |

## Where the cohort comes from

sAMPle's own `sample.benchmark` module reconstructs a fixed, seeded sample of the general
antibacterial classifier's training holdout (`BENCHMARK_SAMPLE_SIZE=1000`, `seed=0`) and compares
every tool at the conventional 0.5 threshold on identical sequences. `external_holdout.fasta` and
`external_holdout_labels.csv` are that exact reconstruction, exported once (from inside the
sAMPle repository, which this package does not depend on at runtime) with:

```python
import sample.benchmark as B

sequences, y_true, entry_ids, df = B._holdout_eval_set(
    "antibacterial/general/any/classification/flaml",
    "data/sample_final_release/final_dataset.csv",
    "data/sample_raw_release/activities.csv",
)
```

## Reproducing `predictions.csv` and the metrics/plots

Only `sample-amp-pred` itself is needed for this part — no sAMPle checkout required:

```bash
pip install .   # from the repository root
sample-amp-pred -i benchmarks/external_holdout.fasta -o /tmp/predictions.csv -f csv -t 0.5
```

Then compare `/tmp/predictions.csv` against `external_holdout_labels.csv` with whatever metrics
library you like (`scikit-learn`'s `roc_auc_score`/`average_precision_score`/etc. is what
`results.json` was computed with). The one sequence with a non-standard residue
(`SAMPLE_DECOY_00049796`, residue `B`) is expected to come back as `prediction=error` rather than
scored — see `results.md` for why that leaves exactly 967 scored sequences.
