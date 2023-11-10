import dendropy
from dendropy.model.discrete import (
    Hky85,
    NucleotideCharacterEvolutionModel,
    simulate_discrete_char_dataset,
    simulate_discrete_chars,
)

from . import pytestmark


def test_qmatrix():
    h = Hky85()

    res = h.qmatrix()
    assert isinstance(res, list)


def test_simulate_discrete_char_dataset():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    m = NucleotideCharacterEvolutionModel()

    res = simulate_discrete_char_dataset(1, t, m)
    assert isinstance(res, dendropy.DataSet)


def test_simulate_discrete_chars():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    m = NucleotideCharacterEvolutionModel()

    res = simulate_discrete_chars(1, t, m)
    assert isinstance(res, dendropy.DnaCharacterMatrix)
