import dendropy
from dendropy.model.birthdeath import (
    birth_death_likelihood,
    birth_death_tree,
    fast_birth_death_tree,
    fit_pure_birth_model,
    fit_pure_birth_model_to_tree,
    uniform_pure_birth_tree,
)

from . import marksmoke as pytestmark


def test_birth_death_tree():
    res = birth_death_tree(1, 1, max_time=1)
    assert isinstance(res, dendropy.Tree)


def test_uniform_pure_birth_tree():
    namespace = dendropy.TaxonNamespace(["A"])

    res = uniform_pure_birth_tree(namespace)
    assert isinstance(res, dendropy.Tree)


def test_fit_pure_birth_model():
    tree = dendropy.Tree()
    tree.seed_node.new_child(edge_length=1)

    res = fit_pure_birth_model(tree=tree)
    assert isinstance(res, dict)


def test_fit_pure_birth_model_to_tree():
    tree = dendropy.Tree()
    tree.seed_node.new_child(edge_length=1)

    res = fit_pure_birth_model_to_tree(tree)
    assert isinstance(res, dict)


def test_birth_death_likelihood():
    tree = dendropy.Tree()
    tree.seed_node.new_child()

    res = birth_death_likelihood(tree=tree, birth_rate=1, death_rate=2)
    assert isinstance(res, float)

def test_fast_birth_death_tree():
    res = fast_birth_death_tree(1, 1, max_time=1)
    assert isinstance(res, dendropy.Tree)
