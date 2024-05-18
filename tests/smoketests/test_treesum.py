import dendropy
from dendropy.calculate.treesum import TopologyCounter, TreeSummarizer
from dendropy.datamodel.treecollectionmodel import SplitDistribution

from . import marksmoke as pytestmark


def test_tree_from_splits():
    dist = SplitDistribution()
    summ = TreeSummarizer()

    res = summ.tree_from_splits(dist)
    assert isinstance(res, dendropy.Tree)


def test_compose_support_label():
    summ = TreeSummarizer()

    res = summ.compose_support_label(1)
    assert isinstance(res, str)


def test_map_split_support_to_node():
    summ = TreeSummarizer()
    tree = dendropy.Tree()
    node = tree.seed_node.new_child()

    res = summ.map_split_support_to_node(node, 1)
    assert isinstance(res, dendropy.Node)


def test_map_split_support_to_tree():
    summ = TreeSummarizer()
    tree = dendropy.Tree()
    dist = SplitDistribution()

    res = summ.map_split_support_to_tree(tree, dist)
    assert isinstance(res, dendropy.Tree)


def test_annotate_nodes_and_edges():
    namespace = dendropy.TaxonNamespace(["A"])
    summ = TreeSummarizer()
    tree = dendropy.Tree(taxon_namespace=namespace)
    dist = SplitDistribution(taxon_namespace=namespace)

    res = summ.annotate_nodes_and_edges(tree, dist)
    assert res is None


def test_summarize_node_ages_on_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A", "B", "C", "D",])
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    tree.seed_node.taxon = taxon_namespace.get_taxon("A")
    tree.seed_node.edge.length = 1
    tree.seed_node.label = "asdf"
    tree.seed_node.age = 1

    tree.update_splits()
    tree.update_bipartitions()

    summ = TreeSummarizer()
    dist = SplitDistribution()

    res = summ.summarize_node_ages_on_tree(tree, dist)
    assert isinstance(res, dendropy.Tree)


def test_summarize_edge_lengths_on_tree():
    tree = dendropy.Tree()
    summ = TreeSummarizer()
    dist = SplitDistribution()

    res = summ.summarize_edge_lengths_on_tree(tree, dist)
    assert isinstance(res, dendropy.Tree)


def test_count_splits_on_trees():
    # TODO: https://github.com/jeetsukumaran/DendroPy/issues/179#issue-1965884280
    namespace = dendropy.TaxonNamespace(["A"])
    summ = TreeSummarizer()
    tree = dendropy.Tree(taxon_namespace=namespace)
    dist = SplitDistribution(taxon_namespace=namespace)
    list = []
    list.append(tree)

    summ.count_splits_on_trees(list, dist)


def test_consensus_tree():
    s = TreeSummarizer()
    t = dendropy.Tree()
    t.seed_node.edge_length = 0.0
    s.consensus_tree([t])


def test_hash_topology():
    c = TopologyCounter()
    t = dendropy.Tree()
    t.update_bipartitions(is_bipartitions_mutable=False)

    res = c.hash_topology(t)
