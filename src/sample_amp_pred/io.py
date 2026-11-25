"""Reading a FASTA file of peptides and writing prediction results in csv/xlsx/json."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import List, Tuple

from .model import PredictionResult

SUPPORTED_FORMATS = ("csv", "xlsx", "json")


def read_fasta(path: str) -> List[Tuple[str, str]]:
    """
    :param path: Path to a FASTA file.
    :returns: ``(identifier, sequence)`` pairs in file order. The identifier is the header line
        up to the first whitespace, with the leading ``>`` removed; unnamed leading sequences
        (a malformed file with no header at all) are named ``peptide_1``, ``peptide_2``, ... .
    :raises ValueError: If the file holds no sequences at all.
    """
    records: List[Tuple[str, str]] = []
    identifier: str | None = None
    chunks: List[str] = []

    def flush() -> None:
        if chunks:
            name = identifier or f"peptide_{len(records) + 1}"
            records.append((name, "".join(chunks)))

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                identifier = line[1:].strip().split()[0] if line[1:].strip() else None
                chunks = []
                continue
            chunks.append(line)
    flush()

    if not records:
        raise ValueError(f"No FASTA sequences found in {path!r}")
    return records


def write_results(results: List[PredictionResult], path: str, fmt: str) -> None:
    """
    :param results: In the order they should be written.
    :param path: Output file path.
    :param fmt: One of :data:`SUPPORTED_FORMATS`.
    """
    fmt = fmt.lower()
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported format {fmt!r}; choose one of {SUPPORTED_FORMATS}")
    if fmt == "csv":
        _write_csv(results, path)
    elif fmt == "json":
        _write_json(results, path)
    else:
        _write_xlsx(results, path)


_COLUMNS = ("id", "sequence", "probability", "prediction", "error")


def _rows(results: List[PredictionResult]) -> List[dict]:
    return [
        {
            "id": r.identifier,
            "sequence": r.sequence,
            "probability": r.probability,
            "prediction": r.prediction,
            "error": r.error or "",
        }
        for r in results
    ]


def _write_csv(results: List[PredictionResult], path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_COLUMNS)
        writer.writeheader()
        writer.writerows(_rows(results))


def _write_json(results: List[PredictionResult], path: str) -> None:
    Path(path).write_text(json.dumps(_rows(results), indent=2), encoding="utf-8")


def _write_xlsx(results: List[PredictionResult], path: str) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Font

    wb = Workbook()
    ws = wb.active
    ws.title = "predictions"
    ws.append(list(_COLUMNS))
    for cell in ws[1]:
        cell.font = Font(bold=True)
    for row in _rows(results):
        ws.append([row[c] for c in _COLUMNS])
    for column_cells in ws.columns:
        length = max(len(str(cell.value)) if cell.value is not None else 0
                     for cell in column_cells)
        ws.column_dimensions[column_cells[0].column_letter].width = min(max(length + 2, 10), 60)
    wb.save(path)
