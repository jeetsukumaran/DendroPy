#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tree test data generation and verification.
"""

import sys
import dendropy

###############################################################################
## Test Tree
# Following tree:
#
#                  a
#                 / \
#                /   \
#               /     \
#              /       \
#             /         \
#            /           \
#           /             c
#          b             / \
#         / \           /   \
#        /   e         /     f
#       /   / \       /     / \
#      /   /   \     g     /   h
#     /   /     \   / \   /   / \
#    i    j     k  l  m  n   o   p
#
#  Can be specified as:
#
#      a -> b -> i;
#      b -> e -> j;
#      e -> k;
#      a -> c;
#      c -> g;
#      c -> f;
#      g -> l;
#      g -> m;
#      f -> n;
#      f -> h -> o;
#      h -> p;
class CuratedTestTree(object):
    dot_str = "a -> b -> i; b -> e -> j; e -> k; a -> c; c -> g; c -> f; g -> l; g -> m; f -> n; f -> h -> o; h -> p;"
    newick_unweighted_edges_str = "((i, (j, k)e)b, ((l, m)g, (n, (o, p)h)f)c)a;"
    newick_weighted_edges_str = "((i:17, (j:3, k:3)e:14)b:33, ((l:6, m:6)g:30, (n:23, (o:11, p:11)h:12)f:13)c:14)a:15;"
    preorder_sequence = ("a", "b", "i", "e", "j", "k", "c", "g", "l", "m", "f", "n", "h", "o", "p",)
    postorder_sequence = ("i", "j", "k", "e", "b", "l", "m", "g", "n", "o", "p", "h", "f", "c", "a",)
    leaf_sequence = ("i", "j", "k", "l", "m", "n", "o", "p",)
    levelorder_sequence = ("a", "b", "c", "i", "e", "g", "f", "j", "k", "l", "m", "n", "h", "o", "p",)
    internal_levelorder_sequence = ("a", "bc", "egf", "h",)
    inorder_sequence = ("i", "b", "j", "e", "k", "a", "l", "g", "m", "c", "n", "f", "o", "h", "p",)
    ageorder_sequence = ("i", "j", "k", "l", "m", "n", "o", "p", "e", "g", "h", "b", "f", "c", "a",)
    leaf_labels = ( "i", "j", "k", "l", "m", "n", "o", "p",)
    internal_labels = ( "a", "b", "c", "e", "f", "g", "h",)
    all_labels = ( "a", "b", "c", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", )
    node_children = {
            "a" : ("b", "c",),
            "b" : ("i", "e",),
            "c" : ("g", "f",),
            "e" : ("j", "k",),
            "f" : ("n", "h",),
            "g" : ("l", "m",),
            "h" : ("o", "p",),
            "i" : (),
            "j" : (),
            "k" : (),
            "l" : (),
            "m" : (),
            "n" : (),
            "o" : (),
            "p" : (),
            }
    node_siblings = {
            "a": (),
            "b": ("c",),
            "c": (),
            "e": (),
            "f": (),
            "g": ("f",),
            "h": (),
            "i": ("e",),
            "j": ("k",),
            "k": (),
            "l": ("m",),
            "m": (),
            "n": ("h",),
            "o": ("p",),
            "p": (),
            }
    node_edge_lengths = {
            "a": 15.0,
            "b": 33.0,
            "c": 14.0,
            "e": 14.0,
            "f": 13.0,
            "g": 30.0,
            "h": 12.0,
            "i": 17.0,
            "j":  3.0,
            "k":  3.0,
            "l":  6.0,
            "m":  6.0,
            "n": 23.0,
            "o": 11.0,
            "p": 11.0,
            }
    node_ages = {
            "a": 50.0,
            "b": 17.0,
            "c": 36.0,
            "e":  3.0,
            "f": 23.0,
            "g":  6.0,
            "h": 11.0,
            "i":  0.0,
            "j":  0.0,
            "k":  0.0,
            "l":  0.0,
            "m":  0.0,
            "n":  0.0,
            "o":  0.0,
            "p":  0.0,
            }
    node_ancestors = {
            "a": (),
            "b": ("a",),
            "c": ("a",),
            "e": ("b", "a",),
            "f": ("c", "a",),
            "g": ("c", "a",),
            "h": ("f", "c", "a",),
            "i": ("b", "a",),
            "j": ("e", "b", "a",),
            "k": ("e", "b", "a",),
            "l": ("g", "c", "a",),
            "m": ("g", "c", "a",),
            "n": ("f", "c", "a",),
            "o": ("h", "f", "c", "a",),
            "p": ("h", "f", "c", "a",),
            }

    # def get_tree(self):
    #     tree = dendropy.Tree()
    #     def add_child_node(parent, label, edge_length):
    #         nd = tree.node_factory()
    #         nd.label = label
    #         nd.edge.length = edge_length
    #         parent.add_child(nd)
    #         return nd
    #     a = tree.seed_node
    #     a.label = "a"
    #     a.edge.length = 15.0
    #     b = add_child_node(a, label="b", edge_length=xxx["b"])
    #     c = add_child_node(a, label="c", edge_length=xxx["c"])
    #     e = add_child_node(b, label="i", edge_length=xxx["i"])
    #     i = add_child_node(b, label="c", edge_length=xxx["c"])
    #     j = add_child_node(e, label="j", edge_length=xxx["j"])
    #     k = add_child_node(e, label="k", edge_length=xxx["k"])
    #     f = add_child_node(c, label="f", edge_length=xxx["f"])
    #     g = add_child_node(c, label="g", edge_length=xxx["g"])
    #     l = add_child_node(g, label="l", edge_length=xxx["l"])
    #     m = add_child_node(g, label="m", edge_length=xxx["m"])
    #     h = add_child_node(f, label="h", edge_length=xxx["h"])
    #     n = add_child_node(f, label="n", edge_length=xxx["n"])
    #     o = add_child_node(h, label="o", edge_length=xxx["o"])
    #     p = add_child_node(h, label="p", edge_length=xxx["p"])
    #     return tree

    def get_tree(self,
            suppress_internal_node_taxa=True,
            suppress_leaf_node_taxa=True,
            taxon_namespace=None,
            node_taxon_label_map=None,
            ):
        if node_taxon_label_map is None:
            node_taxon_label_map = {}
        tree = dendropy.Tree(taxon_namespace=taxon_namespace)
        a = tree.seed_node
        a.label = "a"
        a.edge.length = 15.0
        b = a.new_child(label="b", edge_length=self.node_edge_lengths["b"])
        assert b.label == "b"
        assert b.edge.length == self.node_edge_lengths[b.label]
        assert b.parent_node is a
        assert b.edge.tail_node is a
        assert b in a._child_nodes
        c = a.new_child(label="c", edge_length=self.node_edge_lengths["c"])
        assert c.label == "c"
        assert c.edge.length == self.node_edge_lengths[c.label]
        assert c.parent_node is a
        assert c.edge.tail_node is a
        assert c in a._child_nodes
        i = b.new_child(label="i", edge_length=self.node_edge_lengths["i"])
        assert i.label == "i"
        assert i.edge.length == self.node_edge_lengths[i.label]
        assert i.parent_node is b
        assert i.edge.tail_node is b
        assert i in b._child_nodes
        e = b.new_child(label="e", edge_length=self.node_edge_lengths["e"])
        assert e.label == "e"
        assert e.edge.length == self.node_edge_lengths[e.label]
        assert e.parent_node is b
        assert e.edge.tail_node is b
        assert e in b._child_nodes
        j = e.new_child(label="j", edge_length=self.node_edge_lengths["j"])
        assert j.label == "j"
        assert j.edge.length == self.node_edge_lengths[j.label]
        assert j.parent_node is e
        assert j.edge.tail_node is e
        assert j in e._child_nodes
        k = e.new_child(label="k", edge_length=self.node_edge_lengths["k"])
        assert k.label == "k"
        assert k.edge.length == self.node_edge_lengths[k.label]
        assert k.parent_node is e
        assert k.edge.tail_node is e
        assert k in e._child_nodes
        g = c.new_child(label="g", edge_length=self.node_edge_lengths["g"])
        assert g.label == "g"
        assert g.edge.length == self.node_edge_lengths[g.label]
        assert g.parent_node is c
        assert g.edge.tail_node is c
        assert g in c._child_nodes
        f = c.new_child(label="f", edge_length=self.node_edge_lengths["f"])
        assert f.label == "f"
        assert f.edge.length == self.node_edge_lengths[f.label]
        assert f.parent_node is c
        assert f.edge.tail_node is c
        assert f in c._child_nodes
        l = g.new_child(label="l", edge_length=self.node_edge_lengths["l"])
        assert l.label == "l"
        assert l.edge.length == self.node_edge_lengths[l.label]
        assert l.parent_node is g
        assert l.edge.tail_node is g
        assert l in g._child_nodes
        m = g.new_child(label="m", edge_length=self.node_edge_lengths["m"])
        assert m.label == "m"
        assert m.edge.length == self.node_edge_lengths[m.label]
        assert m.parent_node is g
        assert m.edge.tail_node is g
        assert m in g._child_nodes
        n = f.new_child(label="n", edge_length=self.node_edge_lengths["n"])
        assert n.label == "n"
        assert n.edge.length == self.node_edge_lengths[n.label]
        assert n.parent_node is f
        assert n.edge.tail_node is f
        assert n in f._child_nodes
        h = f.new_child(label="h", edge_length=self.node_edge_lengths["h"])
        assert h.label == "h"
        assert h.edge.length == self.node_edge_lengths[h.label]
        assert h.parent_node is f
        assert h.edge.tail_node is f
        assert h in f._child_nodes
        o = h.new_child(label="o", edge_length=self.node_edge_lengths["o"])
        assert o.label == "o"
        assert o.edge.length == self.node_edge_lengths[o.label]
        assert o.parent_node is h
        assert o.edge.tail_node is h
        assert o in h._child_nodes
        p = h.new_child(label="p", edge_length=self.node_edge_lengths["p"])
        assert p.label == "p"
        assert p.edge.length == self.node_edge_lengths[p.label]
        assert p.parent_node is h
        assert p.edge.tail_node is h
        assert p in h._child_nodes
        tree._debug_check_tree()
        leaf_nodes = set([i, j, k, l, m, n, o, p])
        internal_nodes = set([b, c, e, f, g, h])
        all_nodes = leaf_nodes | internal_nodes | set([a])
        if not suppress_internal_node_taxa:
            for nd in internal_nodes | set([a]):
                label = node_taxon_label_map.get(nd.label, nd.label) # default to same label as node
                t = tree.taxon_namespace.require_taxon(label=label)
                nd.taxon = t
                assert t in tree.taxon_namespace
        if not suppress_leaf_node_taxa:
            for nd in leaf_nodes:
                label = node_taxon_label_map.get(nd.label, nd.label) # default to same label as node
                t = tree.taxon_namespace.require_taxon(label=label)
                nd.taxon = t
                assert t in tree.taxon_namespace
        return tree, all_nodes, leaf_nodes, internal_nodes

    def verify_curated_tree(self,
            tree,
            suppress_internal_node_taxa=True,
            suppress_leaf_node_taxa=False,
            suppress_edge_lengths=False,
            node_taxon_label_map=None):
        if node_taxon_label_map is None:
            node_taxon_label_map = {}
        for nd, exp_nd in zip(tree, self.preorder_sequence):
            if ( (nd.is_leaf() and suppress_leaf_node_taxa)
                    or ((not nd.is_leaf()) and suppress_internal_node_taxa) ):
                label = nd.label
            else:
                self.assertIsNot(nd.taxon, None)
                label = nd.taxon.label
            self.assertEqual(label, node_taxon_label_map.get(exp_nd, exp_nd))
            if suppress_edge_lengths:
                self.assertIs(nd.edge.length, None)
            else:
                self.assertEqual(nd.edge.length, self.node_edge_lengths[exp_nd])
            nd.canonical_label = exp_nd
        for nd in tree:
            children = [c.canonical_label for c in nd.child_node_iter()]
            self.assertCountEqual(children, self.node_children[nd.canonical_label])
            if nd.parent_node is None:
                self.assertEqual(len(self.node_ancestors[nd.canonical_label]), 0)
            else:
                self.assertEqual(nd.parent_node.canonical_label,
                        self.node_ancestors[nd.canonical_label][0])

    def get_newick_string(self,
            suppress_edge_lengths=False,
            node_taxon_label_map=None,
            edge_label_compose_fn=None,
            tree_preamble_tokens=None):
        node_tag = {}
        if node_taxon_label_map is None:
            node_taxon_label_map = {}
        if edge_label_compose_fn is None:
            edge_label_compose_fn = lambda e: "{:6.5E}".format(e)
        node_tag = {}
        for nd in self.preorder_sequence:
            label = node_taxon_label_map.get(nd, nd) # default to same label as node
            if suppress_edge_lengths:
                node_tag[nd] = label
            else:
                node_tag[nd] = "{}:{}".format(label, edge_label_compose_fn(self.node_edge_lengths[nd]))
        if tree_preamble_tokens is None:
            node_tag["preamble"] = ""
        else:
            node_tag["preamble"] = tree_preamble_tokens
        s = "{preamble}(({i}, ({j}, {k}){e}){b}, (({l}, {m}){g}, ({n}, ({o}, {p}){h}){f}){c}){a};".format(
                **node_tag
                )
        return s

