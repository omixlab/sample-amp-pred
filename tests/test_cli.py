import csv
import json

import pytest

from sample_amp_pred.cli import main

FASTA = ">good_one\nGIGKFLHSAKKFGKAFVGEIMNS\n>bad_one\nXJZ\n"


@pytest.fixture
def fasta_file(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(FASTA)
    return path


def test_cli_writes_csv_with_header(tmp_path, fasta_file):
    out = tmp_path / "out.csv"
    assert main(["-i", str(fasta_file), "-o", str(out), "-f", "csv"]) == 0
    rows = list(csv.DictReader(out.open()))
    assert len(rows) == 2
    assert rows[0]["id"] == "good_one"
    assert rows[0]["error"] == ""
    assert float(rows[0]["probability"]) > 0
    assert rows[1]["id"] == "bad_one"
    assert rows[1]["prediction"] == "error"


def test_cli_infers_format_from_extension(tmp_path, fasta_file):
    out = tmp_path / "out.json"
    assert main(["-i", str(fasta_file), "-o", str(out)]) == 0
    data = json.loads(out.read_text())
    assert len(data) == 2


def test_cli_threshold_flag_changes_prediction(tmp_path, fasta_file):
    out_low = tmp_path / "low.csv"
    out_high = tmp_path / "high.csv"
    main(["-i", str(fasta_file), "-o", str(out_low), "-t", "0.0"])
    main(["-i", str(fasta_file), "-o", str(out_high), "-t", "1.0"])
    low_rows = list(csv.DictReader(out_low.open()))
    high_rows = list(csv.DictReader(out_high.open()))
    assert low_rows[0]["prediction"] == "active"
    assert high_rows[0]["prediction"] == "inactive"


def test_cli_missing_input_file_exits_nonzero(tmp_path):
    out = tmp_path / "out.csv"
    assert main(["-i", str(tmp_path / "does_not_exist.fasta"), "-o", str(out)]) != 0


def test_cli_rejects_invalid_threshold(fasta_file, tmp_path, capsys):
    out = tmp_path / "out.csv"
    with pytest.raises(SystemExit):
        main(["-i", str(fasta_file), "-o", str(out), "-t", "1.5"])
