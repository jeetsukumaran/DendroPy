import dendropy
from dendropy.model.treeshape import star_tree

from . import marksmoke as pytestmark


def test_star_tree():
    res = star_tree(taxon_namespace=dendropy.TaxonNamespace(["A"]))
    assert isinstance(res, dendropy.Tree)
