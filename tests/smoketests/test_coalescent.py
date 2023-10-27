from . import pytestmark

from dendropy.model.coalescent import (
    discrete_time_to_coalescence,
    time_to_coalescence,
    expected_tmrca,
    coalesce_nodes,
    node_waiting_time_pairs,
    extract_coalescent_frames,
    log_probability_of_coalescent_frames,
    log_probability_of_coalescent_tree,
    contained_coalescent_tree,
    pure_kingman_tree,
    pure_kingman_tree_shape,
    mean_kingman_tree,
    constrained_kingman_tree,
)

from dendropy.datamodel.treemodel import Tree, Node
import dendropy

def test_discrete_time_to_coalescence():
    # discrete_time_to_coalescence(2) # if no second parameter, errors on 72
    # discrete_time_to_coalescence(1, 1) # first parameter must be >1
    discrete_time_to_coalescence(2, 1)
    assert type(discrete_time_to_coalescence(2, 1, 1)) == int

def test_time_to_coalesence():
    # time_to_coalescence(1) # also requires >1
    time_to_coalescence(2)

def test_expected_tmrca():
    # expected_tmrca(1)
    assert type(expected_tmrca(2)) == float
    assert type(expected_tmrca(2, 1)) == float

def test_coalesce_nodes():
    t = Tree()
    assert type(coalesce_nodes(None) == list)
    assert type(coalesce_nodes(t.nodes(), use_expected_tmrca=True)) == list

def test_node_waiting_time_pairs():
    t = Tree()
    n = Node()
    t.seed_node.add_child(n)
    # requires internal nodes

    assert type(node_waiting_time_pairs(t)) == list

def test_extract_coalescent_frames():
    t = Tree()
    n = Node()
    t.seed_node.add_child(n)
    
    assert type(extract_coalescent_frames(t)) == dict

def test_log_probability_of_coalescent_frames():
    t = Tree()
    t.seed_node.new_child(edge_length=1)
    t.seed_node.new_child(edge_length=1)
    # how did this even pass

    assert type(log_probability_of_coalescent_frames(extract_coalescent_frames(t), 1)) == float

def test_log_probability_of_coalescent_trees():
    t = Tree()
    t.seed_node.new_child(edge_length=1)
    t.seed_node.new_child(edge_length=1)

    log_probability_of_coalescent_tree(t, 1)

def test_contained_coalescent_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = Tree(taxon_namespace=taxon_namespace)

    n = t.seed_node.new_child()
    n.taxon = taxon_namespace.get_taxon("A")

    gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=t.taxon_namespace,
        num_contained=1)

    assert type(contained_coalescent_tree(t, gene_to_species_map)) == Tree

def test_pure_kingman_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    assert type(pure_kingman_tree(taxon_namespace)) == Tree

def test_pure_kingman_tree_shape():
    assert type(pure_kingman_tree_shape(1, 1)) == Tree

def test_mean_kingman_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    assert type(mean_kingman_tree(taxon_namespace)) == Tree

def test_constrained_kingman_tree():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = Tree(taxon_namespace=taxon_namespace)
    n = t.seed_node.new_child()
    n.taxon = taxon_namespace.get_taxon("A")

    assert type(constrained_kingman_tree(t)) == tuple
 