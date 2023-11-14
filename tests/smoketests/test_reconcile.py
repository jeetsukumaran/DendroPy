import dendropy
from dendropy.model.reconcile import ContainingTree

from . import marksmoke as pytestmark


def test_clear():
    # TODO: https://github.com/jeetsukumaran/DendroPy/issues/179
    # taxon_namespace = dendropy.TaxonNamespace(["A"])
    # t = dendropy.Tree(taxon_namespace=taxon_namespace)
    # t.seed_node.taxon = taxon_namespace.get_taxon("A")
    # t.seed_node.label = "A"
    # m = dendropy.TaxonNamespaceMapping(domain_taxon_namespace=taxon_namespace,
    #             range_taxon_namespace=taxon_namespace,
    #             mapping_fn=lambda t: taxon_namespace.require_taxon(label=t.label[0].upper()))
    # c = ContainingTree(t, taxon_namespace, m)

    # res = c.clear()
    # assert res == None
    pass


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
