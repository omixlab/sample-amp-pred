"""
sample-amp-pred -- score peptides for general antibacterial activity with sAMPle's bundled
classifier.

    sample-amp-pred -i peptides.fasta -o predictions.xlsx -f xlsx -t 0.5
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from .io import SUPPORTED_FORMATS, read_fasta, write_results
from .model import DEFAULT_THRESHOLD, load_bundled_model, predict


def _infer_format(output_path: str) -> str | None:
    ext = os.path.splitext(output_path)[1].lstrip(".").lower()
    return ext if ext in SUPPORTED_FORMATS else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sample-amp-pred",
        description="Predict general antibacterial activity of peptides with sAMPle's "
                    "bundled classifier (https://sample-database.org).",
    )
    parser.add_argument("-i", "--input", required=True, metavar="FASTA",
                        help="Input FASTA file of peptide sequences.")
    parser.add_argument("-o", "--output", required=True, metavar="FILE",
                        help="Output file path.")
    parser.add_argument("-f", "--format", choices=SUPPORTED_FORMATS, default=None,
                        help="Output format. Guessed from the output file's extension when "
                            "omitted; defaults to csv if that fails too.")
    parser.add_argument("-t", "--threshold", type=float, default=None, metavar="PROB",
                        help=f"Probability (0-1) at or above which a peptide is reported as "
                            f"'active'. Default: {DEFAULT_THRESHOLD} (sAMPle's own calibrated "
                            f"cutoff for this model).")
    parser.add_argument("--version", action="version", version=f"sample-amp-pred {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.threshold is not None and not (0.0 <= args.threshold <= 1.0):
        parser.error("--threshold must be between 0 and 1")

    fmt = args.format or _infer_format(args.output) or "csv"
    threshold = DEFAULT_THRESHOLD if args.threshold is None else args.threshold

    try:
        records = read_fasta(args.input)
    except (OSError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    print(f"Loaded {len(records)} sequence(s) from {args.input}", file=sys.stderr)

    bundle = load_bundled_model()
    results = predict(records, bundle, threshold)

    n_errors = sum(1 for r in results if r.error)
    n_active = sum(1 for r in results if r.prediction == "active")
    if n_errors:
        print(f"warning: {n_errors} of {len(results)} sequence(s) could not be scored "
              f"(see the 'error' column)", file=sys.stderr)
    print(f"{n_active} of {len(results) - n_errors} scored sequence(s) predicted active "
          f"at threshold {threshold}", file=sys.stderr)

    write_results(results, args.output, fmt)
    print(f"Wrote {args.output} ({fmt})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
