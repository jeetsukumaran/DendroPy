import dendropy

# from . import pytestmark
from dendropy.datamodel.treecollectionmodel import (
    SplitDistribution,
    TreeArray,
    TreeList,
)


def test_maximum_product_of_split_support_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    l = TreeList()
    l.__add__(t)
    l.maximum_product_of_split_support_tree()


def test_maximum_sum_of_split_support_tree():
    t = dendropy.Tree()
    l = TreeList()
    l.__add__(t)
    l.maximum_sum_of_split_support_tree()


def test_frequency_of_split():
    l = TreeList()
    l.frequency_of_split(labels="A")


def test_splits_considered():
    d = SplitDistribution()
    d.splits_considered()


def test_update():
    d = SplitDistribution()
    d1 = SplitDistribution()
    d.update(d1)


def test_collapse_edges_with_less_than_minimum_support():
    d = SplitDistribution()
    t = dendropy.Tree()
    d.collapse_edges_with_less_than_minimum_support(t)


def test_treearray_update():
    a = TreeArray()
    a1 = TreeArray()
    a.update(a1)


def test_extend():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    a = TreeArray(taxon_namespace=taxon_namespace)
    a1 = TreeArray(taxon_namespace=taxon_namespace)
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    a1.add_tree(t)
    a.extend(a1)
