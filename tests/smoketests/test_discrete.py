import dendropy
from dendropy.model.discrete import (
    Hky85,
    NucleotideCharacterEvolutionModel,
    hky85_chars,
    simulate_discrete_char_dataset,
    simulate_discrete_chars,
)

from . import marksmoke as pytestmark


def test_qmatrix():
    h = Hky85()

    res = h.qmatrix()
    assert isinstance(res, list)


def test_simulate_discrete_char_dataset():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    model = NucleotideCharacterEvolutionModel()

    res = simulate_discrete_char_dataset(1, tree, model)
    assert isinstance(res, dendropy.DataSet)


def test_simulate_discrete_chars():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    model = NucleotideCharacterEvolutionModel()

    res = simulate_discrete_chars(1, tree, model)
    assert isinstance(res, dendropy.DnaCharacterMatrix)


def test_hky85_chars():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    hky85_chars(10, tree)
