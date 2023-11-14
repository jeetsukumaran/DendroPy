from dendropy.calculate.combinatorics import (
    choose,
    num_edges_on_tree,
    num_internal_edges_on_tree,
    num_internal_nodes_on_tree,
)

from . import marksmoke as pytestmark


def test_choose():
    res = choose(0, 1)
    assert res == 0


def test_num_edges_on_tree():
    res = num_edges_on_tree(1, False)
    assert isinstance(res, int)


def test_num_internal_edges_on_tree():
    res1 = num_internal_edges_on_tree(1)
    res2 = num_internal_edges_on_tree(1, False)
    assert isinstance(res1, int)
    assert isinstance(res2, int)


def test_num_internal_nodes_on_tree():
    res1 = num_internal_nodes_on_tree(1)
    res2 = num_internal_nodes_on_tree(1, False)
    assert isinstance(res1, int)
    assert isinstance(res2, int)
