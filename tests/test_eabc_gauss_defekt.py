from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_from_lean import EClass, class_of
from collatz_eabc_gauss_defekt_test import (
    classify_prime,
    correlation_report,
    eabc_coarse_bucket,
    gauss_coarse_bucket,
    gauss_split_class,
    run,
)


def test_gauss_split_inert_known_primes():
    assert gauss_split_class(5) == "split"   # 4·1+1
    assert gauss_split_class(13) == "split"  # 4·3+1
    assert gauss_split_class(7) == "inert"   # 4·1+3
    assert gauss_split_class(11) == "inert"  # 4·2+3
    assert gauss_split_class(2) == "ramified"


def test_eabc_class_small_primes():
    assert class_of(5) is EClass.A
    assert class_of(13) is EClass.E
    assert class_of(7) is EClass.B
    assert class_of(11) is EClass.C


def test_coarse_bucket_alignment():
    assert gauss_coarse_bucket("split") == "EA"
    assert gauss_coarse_bucket("inert") == "BC"
    assert eabc_coarse_bucket(EClass.E) == "EA"
    assert eabc_coarse_bucket(EClass.A) == "EA"
    assert eabc_coarse_bucket(EClass.B) == "BC"
    assert eabc_coarse_bucket(EClass.C) == "BC"


def test_classify_prime_row_structure():
    row = classify_prime(13)
    assert row is not None
    assert row.p == 13
    assert row.gauss == "split"
    assert row.eabc == "E"
    assert row.coarse_match is True


def test_correlation_exact_for_small_primes():
    report = correlation_report(200)
    assert report["counts"]["coarse_mismatch"] == 0
    assert report["coarse_match_rate"] == 1.0
    c = report["counts"]
    assert c["split_E"] + c["split_A"] + c["inert_B"] + c["inert_C"] == c["total_p_gt_3"]
    cross = report["cross_table"]
    assert cross["split"]["B"] == 0
    assert cross["split"]["C"] == 0
    assert cross["inert"]["E"] == 0
    assert cross["inert"]["A"] == 0


def test_correlation_exact_up_to_5000():
    report = correlation_report(5000)
    assert report["counts"]["coarse_mismatch"] == 0
    assert report["coarse_match_rate"] == 1.0
    assert "exact_coarse_bipartition" in report["mapping_verdict"]


def test_run_writes_json(tmp_path):
    out = tmp_path / "gauss_defekt.json"
    report = run(max_p=100, output=out)
    assert out.is_file()
    assert report["meta"]["max_p"] == 100
    assert "cross_table" in report
