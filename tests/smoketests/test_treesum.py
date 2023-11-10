import dendropy
from dendropy.calculate.treesum import TreeSummarizer
from dendropy.datamodel.treecollectionmodel import SplitDistribution

from . import pytestmark


def test_tree_from_splits():
    d = SplitDistribution()
    s = TreeSummarizer()

    res = s.tree_from_splits(d)
    assert isinstance(res, dendropy.Tree)


def test_compose_support_label():
    s = TreeSummarizer()

    res = s.compose_support_label(1)
    assert isinstance(res, str)


def test_map_split_support_to_node():
    s = TreeSummarizer()
    t = dendropy.Tree()
    n = t.seed_node.new_child()

    res = s.map_split_support_to_node(n, 1)
    assert isinstance(res, dendropy.Node)


def test_map_split_support_to_tree():
    s = TreeSummarizer()
    t = dendropy.Tree()
    d = SplitDistribution()

    res = s.map_split_support_to_tree(t, d)
    assert isinstance(res, dendropy.Tree)


def test_annotate_nodes_and_edges():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    s = TreeSummarizer()
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    d = SplitDistribution(taxon_namespace=taxon_namespace)

    res = s.annotate_nodes_and_edges(t, d)
    assert res is None


def test_summarize_node_ages_on_tree():
    """
    TODO: https://github.com/jeetsukumaran/DendroPy/issues/179#issue-1965884280

    taxon_namespace = dendropy.TaxonNamespace(["A", "B", "C", "D",])
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    tree.seed_node.taxon = taxon_namespace.get_taxon("A")
    tree.seed_node.edge.length = 1
    tree.seed_node.label = "asdf"
    tree.seed_node.age = 1

    tree.update_splits()
    tree.update_bipartitions()


    s = TreeSummarizer()
    d = SplitDistribution()


    s.summarize_node_ages_on_tree(tree, d)
    """


def test_summarize_edge_lengths_on_tree():
    t = dendropy.Tree()
    s = TreeSummarizer()
    d = SplitDistribution()

    s.summarize_edge_lengths_on_tree(t, d)


def test_count_splits_on_trees():
    """
    TODO: https://github.com/jeetsukumaran/DendroPy/issues/179#issue-1965884280
    s = TreeSummarizer()
    t = dendropy.Tree()
    d = SplitDistribution(taxon_namespace=None)
    l = []
    l.append(t)

    s.count_splits_on_trees(l, d)
    """


def test_consensus_tree():
    """
    TODO: https://github.com/jeetsukumaran/DendroPy/issues/179#issue-1965884280
    s = TreeSummarizer()
    t = dendropy.Tree()
    l = []
    l.append(t)

    s.consensus_tree(l)
    """


def test_hash_topology():
    """AssertionError: Bipartition is mutable: hash is unstable
    c = TopologyCounter()
    t = dendropy.Tree()
    t.update_bipartitions(is_bipartitions_mutable=False)

    res = c.hash_topology(t)
    """
