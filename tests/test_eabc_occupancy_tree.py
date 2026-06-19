"""Tests fuer EABC_OCCUPANCY_TREE (merge_state, C_M, Referenz N=10^8)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eabc_occupancy_tree import (
    M,
    admissible_quad_channels,
    channel_universe,
    eabc_fiber,
    identity_state,
    merge_state,
    occupancy_tree,
    reduce_tree,
    scan_occupancy,
    state_equal,
)


def test_channel_universe_size():
    assert len(channel_universe()) == 378
    assert len(admissible_quad_channels()) == 378


def test_occupancy_tree_api_small():
    result = occupancy_tree(10_000, block_size=5_000)
    assert result["channels_total"] == 378
    assert 0 < result["occupied"] <= 378
    assert 0.0 < result["rho"] <= 1.0
    assert set(result["fiber_total"]) == {"A", "C"}


def test_identity_state_is_neutral():
    Z0 = identity_state()
    Z1 = ({1, 2}, {1: 10, 2: 20}, {1: 3, 2: 1})
    assert state_equal(merge_state(Z0, Z1), Z1)
    assert state_equal(merge_state(Z1, Z0), Z1)


def test_merge_state_associative():
    Z1 = ({1, 2}, {1: 10, 2: 20}, {1: 3, 2: 1})
    Z2 = ({2, 3}, {2: 15, 3: 30}, {2: 2, 3: 5})
    Z3 = ({4}, {4: 40}, {4: 7})
    left = merge_state(merge_state(Z1, Z2), Z3)
    right = merge_state(Z1, merge_state(Z2, Z3))
    assert state_equal(left, right)
    O, T, n = left
    assert O == {1, 2, 3, 4}
    assert T == {1: 10, 2: 15, 3: 30, 4: 40}
    assert n == {1: 3, 2: 3, 3: 5, 4: 7}


def test_merge_state_commutative():
    Z1 = ({1, 2}, {1: 10, 2: 20}, {1: 3, 2: 1})
    Z2 = ({2, 3}, {2: 15, 3: 30}, {2: 2, 3: 5})
    assert state_equal(merge_state(Z1, Z2), merge_state(Z2, Z1))


def test_reduce_tree_matches_iterative_merge():
    states = [
        ({1}, {1: 1}, {1: 1}),
        ({2}, {2: 2}, {2: 1}),
        ({3}, {3: 3}, {3: 1}),
    ]
    tree = reduce_tree(states)
    manual = states[0]
    for s in states[1:]:
        manual = merge_state(manual, s)
    assert tree == manual


def test_eabc_fiber_labels():
    assert eabc_fiber(1) == "E"
    assert eabc_fiber(5) == "A"
    assert eabc_fiber(7) == "B"
    assert eabc_fiber(11) == "C"
    assert eabc_fiber(50381) == "A"


def test_scan_occupancy_reference_1e8():
    report = scan_occupancy(10**8, block_width=2_000_000)
    assert report["card_C_M"] == 378
    assert report["card_O_N"] == 378
    assert abs(report["rho_N"] - 1.0) < 1e-12
    assert report["total_events"] == 4766
    assert report["latest_T"]["T_max"] == 54_044_321
    assert report["latest_T"]["channels"] == [50_381]
    assert report["latest_T"]["fiber"] == "A"