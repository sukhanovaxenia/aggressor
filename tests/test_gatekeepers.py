"""Unit tests for gatekeeper slot detection and ranking."""
from aggressor.analysis.gatekeepers import (
    select_gatekeeper_positions,
    best_gatekeeper_for_position,
    rank_gatekeeper_slots,
    flanking_gatekeeper_positions,
)

GK = {"P", "K", "R", "D", "E"}


def test_auto_mode_returns_both_boundaries():
    assert select_gatekeeper_positions([5, 6, 7], "AAAAIFLAA") == {5, 7}


def test_cap_keeps_highest_propensity_slot():
    # I (pos 5, propensity 1.822) outranks L (pos 7, 1.380)
    assert select_gatekeeper_positions([5, 6, 7], "AAAAIFLAA", max_gatekeepers=1) == {5}


def test_best_gatekeeper_is_lowest_propensity():
    assert best_gatekeeper_for_position(5, "AAAAIFLAA", GK) == "K"  # K most suppressive


def test_already_gatekeeper_returns_none():
    assert best_gatekeeper_for_position(1, "KAAAA", GK) is None


def test_ranking_is_descending_by_gain():
    ranked = rank_gatekeeper_slots([5, 7], "AAAAIFLAA", GK)
    gains = [g for _, g in ranked]
    assert gains == sorted(gains, reverse=True)


def test_flanks_are_clamped_to_region():
    # Cluster at the region edge: left flank clipped, right flank kept
    assert flanking_gatekeeper_positions([1, 2, 3], region_start=1, region_end=10) == [4]
