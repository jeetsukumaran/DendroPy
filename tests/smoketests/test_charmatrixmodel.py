from types import GeneratorType

import dendropy
from dendropy.datamodel.charmatrixmodel import CharacterMatrix, DiscreteCharacterMatrix

from . import marksmoke as pytestmark


def test_concatenate():
    namespace = dendropy.TaxonNamespace(["A"])
    mat = CharacterMatrix()
    mat.taxon_namespace = namespace
    mat["A"]

    res = CharacterMatrix.concatenate([mat])
    assert isinstance(res, CharacterMatrix)


def test_vectors():
    mat = CharacterMatrix()

    res = mat.vectors()
    assert isinstance(res, list)


def test_items():
    mat = CharacterMatrix()

    res = mat.items()
    assert isinstance(res, GeneratorType)


def test_remove_sequences():
    namespace = dendropy.TaxonNamespace(["A"])
    mat = CharacterMatrix()
    mat.taxon_namespace = namespace

    res = mat.remove_sequences([])
    assert res is None


def test_discard_sequences():
    mat = CharacterMatrix()

    res = mat.discard_sequences([])
    assert res is None


def test_keep_sequences():
    mat = CharacterMatrix()

    res = mat.keep_sequences([])
    assert res is None


def test_description():
    mat = CharacterMatrix()

    res = mat.description()
    assert isinstance(res, str)


def test_append_taxon_sequence():
    namespace = dendropy.TaxonNamespace(["A"])
    mat = DiscreteCharacterMatrix(taxon_namespace=namespace)
    res = mat.append_taxon_sequence("A", [])
    assert res is None


def test_remap_to_state_alphabet_by_symbol():
    mat = DiscreteCharacterMatrix()

    res = mat.remap_to_state_alphabet_by_symbol("A")
    assert res is None
