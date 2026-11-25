"""
Loads the bundled sAMPle general antibacterial classifier and scores peptide sequences with it.

The model file (``data/general_antibacterial_model.pkl``) is a plain dict -- a fitted
``xgboost.XGBClassifier``, its training column order, and its calibrated decision threshold --
exported from sAMPle's own served model artifact. It carries no dependency on the sAMPle
codebase itself, only on ``xgboost``/``scikit-learn`` (to unpickle the classifier) and ``numpy``.
"""
from __future__ import annotations

import importlib.resources
import pickle
from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

from . import features

MODEL_RESOURCE = ("sample_amp_pred.data", "general_antibacterial_model.pkl")

#: sAMPle's own calibrated operating point for this model, used when the caller passes no
#: threshold of their own. This is the general antibacterial classifier specifically: it predates
#: the out-of-fold threshold calibration sAMPle's pathogen-specific classifiers use (see that
#: project's own manuscript), so its trained threshold is the conventional 0.5 -- not a value
#: this tool invented.
DEFAULT_THRESHOLD = 0.5


@dataclass(frozen=True)
class ModelBundle:
    estimator: Any
    columns: Tuple[str, ...]
    threshold: float
    namespace: str
    holdout_metrics: Dict[str, float]


def load_bundled_model() -> ModelBundle:
    """
    :returns: The bundled model, ready to score feature matrices built by
        :func:`sample_amp_pred.features.build_feature_matrix` with ``.columns``.
    :raises FileNotFoundError: If the package was installed without its data file (should not
        happen through a normal ``pip install`` of this distribution).
    """
    package, filename = MODEL_RESOURCE
    data = importlib.resources.files(package).joinpath(filename).read_bytes()
    bundle = pickle.loads(data)
    return ModelBundle(
        estimator=bundle["estimator"],
        columns=tuple(bundle["columns"]),
        threshold=float(bundle.get("threshold", DEFAULT_THRESHOLD)),
        namespace=str(bundle.get("namespace", "antibacterial/general/any/classification")),
        holdout_metrics=dict(bundle.get("holdout_metrics", {})),
    )


@dataclass(frozen=True)
class PredictionResult:
    identifier: str
    sequence: str
    probability: float | None
    prediction: str
    error: str | None = None


def predict(records: Sequence[Tuple[str, str]], bundle: ModelBundle,
           threshold: float) -> List[PredictionResult]:
    """
    Scores `records` against the general antibacterial classifier.

    :param records: ``(identifier, raw_sequence)`` pairs, in the order they should be reported.
    :param bundle: From :func:`load_bundled_model`.
    :param threshold: Probability at or above which a peptide is reported as ``active``.
    :returns: One :class:`PredictionResult` per record, in the same order. A record whose
        sequence cannot be scored (empty, too short, or carrying a non-standard residue) gets
        ``probability=None`` and its ``error`` set, rather than being silently dropped -- every
        input FASTA entry gets exactly one output row.
    """
    scorable: List[int] = []
    cleaned: List[str] = []
    results: List[PredictionResult] = [None] * len(records)  # type: ignore[list-item]

    for i, (identifier, raw_sequence) in enumerate(records):
        seq = features.clean_sequence(raw_sequence)
        problem = features.validate_sequence(seq)
        if problem is not None:
            results[i] = PredictionResult(identifier, raw_sequence, None, "error", problem)
            continue
        scorable.append(i)
        cleaned.append(seq)

    if cleaned:
        matrix = features.build_feature_matrix(cleaned, bundle.columns)
        probabilities = np.asarray(bundle.estimator.predict_proba(matrix))[:, 1]
        for pos, i in enumerate(scorable):
            probability = float(probabilities[pos])
            label = "active" if probability >= threshold else "inactive"
            results[i] = PredictionResult(records[i][0], records[i][1], probability, label)

    return results  # type: ignore[return-value]
