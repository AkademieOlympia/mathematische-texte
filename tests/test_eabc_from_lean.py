from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_from_lean import (
    GAP_SIGNATURE,
    Chirality,
    EClass,
    center,
    chirality_order,
    class_of,
    definitions_overview,
    is_prime_quadruplet,
    q,
    residue,
    t,
    t4,
)


def test_residue_mapping_matches_lean():
    assert residue(EClass.E) == 1
    assert residue(EClass.A) == 5
    assert residue(EClass.B) == 7
    assert residue(EClass.C) == 11


def test_class_of_for_representative_values():
    assert class_of(1) is EClass.E
    assert class_of(5) is EClass.A
    assert class_of(7) is EClass.B
    assert class_of(11) is EClass.C
    assert class_of(13) is EClass.E
    assert class_of(17) is EClass.A
    assert class_of(0) is None
    assert class_of(2) is None


def test_quadruplet_structure_center_and_gaps():
    values = q(5)
    assert values == [5, 7, 11, 13]
    assert center(5) == 9
    assert GAP_SIGNATURE == [2, 4, 2]
    assert [values[i + 1] - values[i] for i in range(3)] == GAP_SIGNATURE


def test_is_prime_quadruplet_examples():
    assert is_prime_quadruplet(5)
    assert not is_prime_quadruplet(7)


def test_chirality_orders():
    assert chirality_order(Chirality.ABCE) == [EClass.A, EClass.B, EClass.C, EClass.E]
    assert chirality_order(Chirality.CEAB) == [EClass.C, EClass.E, EClass.A, EClass.B]


def test_rotation_operator_t_and_t4():
    assert t(EClass.E) is EClass.A
    assert t(EClass.A) is EClass.B
    assert t(EClass.B) is EClass.C
    assert t(EClass.C) is EClass.E

    for x in EClass:
        assert t4(x) is x


def test_class_of_mod12_property_on_large_range():
    expected = {
        1: EClass.E,
        5: EClass.A,
        7: EClass.B,
        11: EClass.C,
    }
    for n in range(0, 10_000):
        assert class_of(n) is expected.get(n % 12)


def test_residue_roundtrip_property():
    for cls in EClass:
        assert class_of(residue(cls)) is cls
        assert class_of(residue(cls) + 12 * 100) is cls


def test_rotation_cycle_lengths():
    for x in EClass:
        y1 = t(x)
        y2 = t(y1)
        y3 = t(y2)
        y4 = t(y3)
        assert y4 is x
        assert y1 is not x
        assert y2 is not x
        assert y3 is not x


def _is_prime_reference(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def test_random_prime_quadruplet_consistency_with_q():
    rng = random.Random(20260527)
    for _ in range(300):
        p = rng.randint(2, 5000)
        vals = q(p)
        assert vals == [p, p + 2, p + 6, p + 8]
        expected = all(_is_prime_reference(v) for v in vals)
        assert is_prime_quadruplet(p) == expected


def test_random_center_property():
    rng = random.Random(314159)
    for _ in range(300):
        p = rng.randint(0, 100_000)
        vals = q(p)
        assert center(p) * 2 == vals[0] + vals[3]


def test_definitions_overview_structure():
    overview = definitions_overview()
    assert set(overview.keys()) == {"enums", "constants", "functions"}
    assert overview["enums"]["EClass"]["values"] == ["E", "A", "B", "C"]
    assert overview["enums"]["Chirality"]["values"] == ["ABCE", "CEAB"]
    assert overview["constants"]["GAP_SIGNATURE"]["value"] == [2, 4, 2]
    assert "class_of" in overview["functions"]
