from sample_amp_pred import features


def test_column_counts():
    assert len(features.DESCRIPTOR_COLUMNS) == 11
    assert len(features.COMPOSITION_COLUMNS) == 28
    assert len(features.DPC_COLUMNS) == 400
    assert len(features.CTD_COLUMNS) == 147
    assert len(features.ALL_COLUMNS) == 586
    assert len(set(features.ALL_COLUMNS)) == 586  # no duplicate column names


def test_clean_sequence_strips_non_letters():
    assert features.clean_sequence(" gigk-flh\n2akk ") == "GIGKFLHAKK"


def test_validate_sequence_rejects_empty():
    assert features.validate_sequence("") is not None


def test_validate_sequence_rejects_too_short():
    assert features.validate_sequence("A") is not None
    assert features.validate_sequence("AA") is None


def test_validate_sequence_rejects_non_standard_residues():
    problem = features.validate_sequence("GIGKXJZ")
    assert problem is not None
    assert "X" in problem and "J" in problem and "Z" in problem


def test_build_feature_matrix_shape_and_order():
    matrix = features.build_feature_matrix(["GIGKFLHSAKKFGKAFVGEIMNS"], features.ALL_COLUMNS)
    assert matrix.shape == (1, 586)


def test_composition_row_sums_to_one():
    row = features._composition_row("GIGKFLHSAKKFGKAFVGEIMNS")
    # the first 20 values are single-residue frequencies, which must sum to 1
    assert abs(sum(row[:20]) - 1.0) < 1e-9


def test_dpc_row_sums_to_one():
    row = features._dpc_row("GIGKFLHSAKKFGKAFVGEIMNS")
    assert abs(sum(row) - 1.0) < 1e-9


def test_dpc_row_all_zero_below_two_residues():
    assert features._dpc_row("A") == [0.0] * len(features.DPC_COLUMNS)
