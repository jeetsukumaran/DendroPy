import tempfile

import dendropy
from dendropy.datamodel.taxonmodel import Taxon, TaxonNamespace, TaxonNamespacePartition
from dendropy.datamodel.treemodel import Bipartition

from . import marksmoke as pytestmark


def test_taxa_bipartition():
    namespace = TaxonNamespace(["A"])

    res = namespace.taxa_bipartition(labels=[])
    assert isinstance(res, Bipartition)


def test_bitmask_taxa_list():
    namespace = TaxonNamespace()

    res = namespace.bitmask_taxa_list([])
    assert isinstance(res, list)


def test_description():
    namespace = TaxonNamespace()

    res = namespace.description()
    assert isinstance(res, str)


def test_taxon_description():
    taxon = Taxon()

    res = taxon.description()
    assert isinstance(res, str)


def test_apply_membership_lists():
    namespace = TaxonNamespace(["A"])
    part = TaxonNamespacePartition(namespace)

    res = part.apply_membership_lists([])
    assert isinstance(res, set)


def test_apply_mapping_attr_name():
    namespace = TaxonNamespace(["A"])
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )

    res = map.apply_mapping_attr_name("label", "A")
    assert res is None


def test_apply_mapping_dict():
    namespace = TaxonNamespace(["A"])
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )

    res = map.apply_mapping_dict({})
    assert res is None


def test_mesquite_association_rows():
    namespace = TaxonNamespace(["A"])
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )

    res = map.mesquite_association_rows()
    assert isinstance(res, str)


def test_write_mesquite_association_block():
    namespace = TaxonNamespace(["A"])
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )

    with tempfile.NamedTemporaryFile(mode="w") as file:
        res = map.write_mesquite_association_block(file)
    assert res is None
