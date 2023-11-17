import dendropy
from dendropy.model.reconcile import ContainingTree

from . import marksmoke as pytestmark


def test_clear():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    tree.seed_node.label = "A"
    map = dendropy.TaxonNamespaceMapping(domain_taxon_namespace=namespace,
                range_taxon_namespace=namespace,
                mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()))
    cont = ContainingTree(tree, namespace, map)

    res = cont.clear()
    assert res is None


def test_rebuild():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )
    cont = ContainingTree(tree, namespace, map)

    res = cont.rebuild()
    assert res is None


def test_embed_contained_kingman():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )
    cont = ContainingTree(tree, namespace, map)

    res = cont.embed_contained_kingman()
    assert isinstance(res, dendropy.Tree)


def test_simulate_contained_kingman():
    namespace = dendropy.TaxonNamespace(["A"])
    tree = dendropy.Tree(taxon_namespace=namespace)
    tree.seed_node.taxon = namespace.get_taxon("A")
    map = dendropy.TaxonNamespaceMapping(
        domain_taxon_namespace=namespace,
        range_taxon_namespace=namespace,
        mapping_fn=lambda t: namespace.require_taxon(label=t.label[0].upper()),
    )
    cont = ContainingTree(tree, namespace, map)

    res = cont.simulate_contained_kingman()
    assert isinstance(res, dendropy.Tree)
