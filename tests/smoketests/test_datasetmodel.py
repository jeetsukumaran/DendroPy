from dendropy.datamodel.datasetmodel import DataSet
import dendropy

from . import marksmoke as pytestmark


def test_unify_taxon_namespaces():
    namespace = dendropy.TaxonNamespace(["A"])
    set = DataSet()
    
    res = set.unify_taxon_namespaces(namespace)
    assert res is None
    
