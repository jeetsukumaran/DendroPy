import dendropy
from dendropy.model.reconcile import ContainingTree

from . import pytestmark


def test_clear():
    """TODO: https://github.com/jeetsukumaran/DendroPy/issues/179
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    t.seed_node.label = "A"
    m = dendropy.TaxonNamespaceMapping(domain_taxon_namespace=taxon_namespace,
                range_taxon_namespace=taxon_namespace,
                mapping_fn=lambda t: taxon_namespace.require_taxon(label=t.label[0].upper()))
    c = ContainingTree(t, taxon_namespace, m)

    res = c.clear()
    assert res == None
    """


def test_rebuild():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    m = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=taxon_namespace,
        range_taxon_namespace=taxon_namespace,
        mapping_fn=lambda t: taxon_namespace.require_taxon(label=t.label[0].upper()),
    )
    c = ContainingTree(t, taxon_namespace, m)

    res = c.rebuild()
    assert res is None


def test_embed_contained_kingman():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    m = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=taxon_namespace,
        range_taxon_namespace=taxon_namespace,
        mapping_fn=lambda t: taxon_namespace.require_taxon(label=t.label[0].upper()),
    )
    c = ContainingTree(t, taxon_namespace, m)

    res = c.embed_contained_kingman()
    assert isinstance(res, dendropy.Tree)


def test_simulate_contained_kingman():
    taxon_namespace = dendropy.TaxonNamespace(["A"])
    t = dendropy.Tree(taxon_namespace=taxon_namespace)
    t.seed_node.taxon = taxon_namespace.get_taxon("A")
    m = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=taxon_namespace,
        range_taxon_namespace=taxon_namespace,
        mapping_fn=lambda t: taxon_namespace.require_taxon(label=t.label[0].upper()),
    )
    c = ContainingTree(t, taxon_namespace, m)

    res = c.simulate_contained_kingman()
    assert isinstance(res, dendropy.Tree)
