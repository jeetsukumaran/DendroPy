import dendropy
from dendropy.datamodel.treecollectionmodel import (
    SplitDistribution,
    TreeArray,
    TreeList,
)

from . import marksmoke as pytestmark


def test_maximum_product_of_split_support_tree():
    list = TreeList.get(data="((A,B),(C,D));((A,C),(B,D));", schema="newick")

    res = list.maximum_product_of_split_support_tree()
    assert isinstance(res, dendropy.Tree)


def test_maximum_sum_of_split_support_tree():
    list = TreeList.get(data="((A,B),(C,D));((A,C),(B,D));", schema="newick")

    res = list.maximum_sum_of_split_support_tree()
    assert isinstance(res, dendropy.Tree)


def test_frequency_of_split():
    list = TreeList.get(data="((A,B),(C,D));((A,C),(B,D));", schema="newick")

    res = list.frequency_of_split(labels="A")
    assert isinstance(res, float)


def test_splits_considered():
    dist = SplitDistribution()

    res = dist.splits_considered()
    assert isinstance(res, tuple)


def test_update():
    dist1 = SplitDistribution()
    dist2 = SplitDistribution()

    res = dist1.update(dist2)
    assert res is None


def test_collapse_edges_with_less_than_minimum_support():
    dist = SplitDistribution()
    tree = dendropy.Tree()

    res = dist.collapse_edges_with_less_than_minimum_support(tree)
    assert res is None


def test_treearray_update():
    arr1 = TreeArray()
    arr2 = TreeArray()

    res = arr1.update(arr2)
    assert res is None


def test_extend():
    namespace = dendropy.TaxonNamespace(["A"])
    arr1 = TreeArray(taxon_namespace=namespace)
    arr2 = TreeArray(taxon_namespace=namespace)

    res = arr1.extend(arr2)
    assert isinstance(res, TreeArray)
