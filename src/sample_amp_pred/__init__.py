"""sample-amp-pred: predict peptide antibacterial activity with sAMPle's bundled classifier."""

__version__ = "0.1.0"

from .model import DEFAULT_THRESHOLD, ModelBundle, PredictionResult, load_bundled_model, predict

__all__ = [
    "__version__",
    "DEFAULT_THRESHOLD",
    "ModelBundle",
    "PredictionResult",
    "load_bundled_model",
    "predict",
]
