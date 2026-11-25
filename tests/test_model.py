from sample_amp_pred.model import DEFAULT_THRESHOLD, load_bundled_model, predict

MAGAININ_2 = "GIGKFLHSAKKFGKAFVGEIMNS"


def test_load_bundled_model():
    bundle = load_bundled_model()
    assert len(bundle.columns) == 586
    assert 0.0 <= bundle.threshold <= 1.0


def test_predict_known_amp_scores_high():
    bundle = load_bundled_model()
    results = predict([("magainin_2", MAGAININ_2)], bundle, DEFAULT_THRESHOLD)
    assert len(results) == 1
    r = results[0]
    assert r.error is None
    assert r.probability is not None
    assert 0.0 <= r.probability <= 1.0
    assert r.prediction in ("active", "inactive")


def test_predict_reports_one_result_per_input_in_order():
    bundle = load_bundled_model()
    records = [("a", MAGAININ_2), ("b", "XJZ"), ("c", "AA")]
    results = predict(records, bundle, DEFAULT_THRESHOLD)
    assert [r.identifier for r in results] == ["a", "b", "c"]
    assert results[0].error is None
    assert results[1].error is not None  # non-standard residues
    assert results[2].error is None      # "AA" is a valid 2-residue sequence


def test_threshold_changes_label_not_probability():
    bundle = load_bundled_model()
    low = predict([("x", MAGAININ_2)], bundle, 0.0)[0]
    high = predict([("x", MAGAININ_2)], bundle, 1.0)[0]
    assert low.probability == high.probability
    assert low.prediction == "active"
    assert high.prediction == "inactive"
