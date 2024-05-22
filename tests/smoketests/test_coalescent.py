import dendropy
from dendropy.model.coalescent import (
    coalesce_nodes,
    constrained_kingman_tree,
    contained_coalescent_tree,
    discrete_time_to_coalescence,
    expected_tmrca,
    extract_coalescent_frames,
    log_probability_of_coalescent_frames,
    log_probability_of_coalescent_tree,
    mean_kingman_tree,
    node_waiting_time_pairs,
    pure_kingman_tree,
    pure_kingman_tree_shape,
    time_to_coalescence,
)

from . import marksmoke as pytestmark


def test_discrete_time_to_coalescence():
    res1 = discrete_time_to_coalescence(2, 1)
    res2 = discrete_time_to_coalescence(2, 1, 1)
    assert isinstance(res1, int)
    assert isinstance(res2, int)


def test_time_to_coalesence():
    res = time_to_coalescence(2)
    assert isinstance(res, float)


def test_expected_tmrca():
    res1 = expected_tmrca(2)
    res2 = expected_tmrca(2, 1)
    assert isinstance(res1, float)
    assert isinstance(res2, float)


def test_coalesce_nodes():
    tree = dendropy.Tree()

    res1 = coalesce_nodes(None)
    res2 = coalesce_nodes(tree.nodes(), use_expected_tmrca=True)
    assert isinstance(res1, list)
    assert isinstance(res2, list)


def test_node_waiting_time_pairs():
    tree = dendropy.Tree()
    node = dendropy.Node()
    tree.seed_node.add_child(node)

    res = node_waiting_time_pairs(tree)
    assert isinstance(res, list)


def test_extract_coalescent_frames():
    tree = dendropy.Tree()
    node = dendropy.Node()
    tree.seed_node.add_child(node)

    res = extract_coalescent_frames(tree)
    assert isinstance(res, dict)


def test_log_probability_of_coalescent_frames():
    tree = dendropy.Tree()
    tree.seed_node.new_child(edge_length=1)
    tree.seed_node.new_child(edge_length=1)

    res = log_probability_of_coalescent_frames(extract_coalescent_frames(tree), 1)
    assert isinstance(res, float)


def test_log_probability_of_coalescent_trees():
    tree = dendropy.Tree()
    tree.seed_node.new_child(edge_length=1)
    tree.seed_node.new_child(edge_length=1)

    res = log_probability_of_coalescent_tree(tree, 1)
    assert isinstance(res, float)


def test_contained_coalescent_tree():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)

    node = tree.seed_node.new_child()
    node.taxon = namespace.get_taxon("A")

    gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=tree.taxon_namespace, num_contained=1
    )

    res = contained_coalescent_tree(tree, gene_to_species_map)
    assert isinstance(res, dendropy.Tree)


def test_pure_kingman_tree():
    namespace = dendropy.TaxonNamespace(["A"])

    res = pure_kingman_tree(namespace)
    assert isinstance(res, dendropy.Tree)


def test_pure_kingman_tree_shape():
    res = pure_kingman_tree_shape(1, 1)
    assert isinstance(res, dendropy.Tree)


def test_mean_kingman_tree():
    namespace = dendropy.TaxonNamespace(["A"])

    res = mean_kingman_tree(namespace)
    assert isinstance(res, dendropy.Tree)


def test_constrained_kingman_tree():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    node = tree.seed_node.new_child()
    node.taxon = namespace.get_taxon("A")

    res = constrained_kingman_tree(tree)
    assert isinstance(res, tuple)
