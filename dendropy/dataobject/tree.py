#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from cStringIO import StringIO
import copy
import re

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

from dendropy.utility import iosys
from dendropy.utility import errors
from dendropy.utility import texttools
from dendropy.dataobject.base import IdTagged
from dendropy.dataobject.taxon import TaxonSetLinked, TaxonLinked

##############################################################################
## TreeList

class TreeList(list, TaxonSetLinked, iosys.Readable, iosys.Writeable):
    """
    Collects and coordinates a list of trees with the associated with the
    same set of taxa.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes a new TreeList object, populating it with any iterable
        container with Tree object members passed as unnamed argument, or
        from a data source if `stream` and `format` are passed.

        TreeList objects can thus be instantiated in the following ways::

            # /usr/bin/env python

            import StringIO
            import dendropy

            # empty tree
            tlst1 = TreeList()

            # populated from list of Tree objects
            t1 = Tree(stream=StringIO("((A,B),(C,D))"), format="newick")
            t2 = Tree(stream=StringIO("((A,C),(B,D))"), format="newick")
            tlist2 = TreeList([t1, t2])

            # tree from data source
            tlst3 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"), format="newick") # same

            # passing keywords to underlying tree parser
            tlst4 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"),
                             format="newick",
                             taxon_set=tlst3.taxon_set,
                             encode_splits=True)

            # shallow-copied from another tree list
            tlst5 = TreeList(t4)

            # deep-copied (but shallow-copy taxa) from another tree list
            tls6 = TreeList([Tree(t) for t in tlst5])

            # can also call `read()` on a TreeList object
            tlst7 = TreeList()
            tlst7.read(StringIO("((A,B),(C,D));((A,C),(B,D));"), "newick")
            tlst7.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")
            tlst7.read_from_path("boot.tre", "newick")

        """
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        iosys.Readable.__init__(self)
        iosys.Writeable.__init__(self)
        list.__init__(self)
        if len(args) > 1:
            raise errors.TooManyArgumentsError(self.__class__.__name__, 1, args)
        elif len(args) == 1:
            if "stream" in kwargs or "format" in kwargs:
                raise errors.MultipleInitializationSourceError(self.__class__.__name__, args[0])
            if hasattr(args[0], "__iter__") and not isinstance(args[0], str):
                for t in args[0]:
                    if not isinstance(t, Tree):
                        raiseTypeError("TreeList() only accepts Tree objects as members")
                    self.append(t)
            else:
                raise errors.InvalidArgumentTypeError(self.__class__.__name__, args[0])
        else:
            self.process_source_kwargs(**kwargs)

        if "oid" in kwargs:
            self.oid = kwargs["oid"]
        if "label" in kwargs:
            self.label = kwargs["label"]

    def __deepcopy__(self, memo):
        # we treat the taxa as immutable and copy the reference even in a deepcopy
        o = self.__class__(taxon_set=self.taxon_set, label=self.label)
        memo[id(self)] = o
        memo[id(self.taxon_set)] = o.taxon_set
        for i, t in enumerate(self.taxon_set):
            memo[id(t)] = o.taxon_set[i]
        for k, v in self.__dict__.iteritems():
            if k not in ['taxon_set', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

    def __str__(self):
        return "[%s]" % " ".join([(str(t)+";") for t in self])

    def read(self, stream, format, **kwargs):
        """
        Populates the `TreeList` from a `format`-formatted file-like
        source `stream`. `format` must be a recognized and tree file
        format, such as `nexus`, `newick`, etc, for which a specialized
        tree list writer is available. If this is not implemented for
        the format specified, then a `UnsupportedFormatError` is raised.

        The following optional keyword arguments are recognized:

            - `encode_splits` specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `translate_dict` should provide a dictionary mapping taxon numbers
               (as found in the source) to taxon labels (as defined in the source).
            - `rooted` specifies the default rooting interpretation of the tree
               (see `dendropy.dataio.nexustokenizer` for details).
            - `finish_node_func` is a function that will be applied to each node
               after it has been constructed.
            - `edge_len_type` specifies the type of the edge lengths (int or float)

        Other keyword arguments may be available, depending on the implementation
        of the reader specialized to handle `format` formats.
        """
        from dendropy.utility import iosys
        from dendropy.dataio import tree_source_iter
        if "taxon_set" in kwargs:
            if kwargs["taxon_set"] is not self.taxon_set:
                raise Exception("Cannot specify a different TaxonSet when reading into an existing TreeList.")
        else:
            kwargs["taxon_set"] = self.taxon_set
        for t in tree_source_iter(stream=stream, format=format, **kwargs):
            if t is not None:
                self.append(t)

    def write(self, **kwargs):
        """
        Writes out `TreeList` in `format` to a destination described by
        one of: `file` or `path`:

            - `file`: A file- or file-like object.
            - `path`: A string specifying the path to a file.

        `format` must be a recognized and tree file format, such as `nexus`,
        `newick`, etc, for which a specialized tree list writer is
        available. If this is not implemented for the format specified, then
        a `UnsupportedFormatError` is raised.

        Additionally, for some formats, the following keywords are recognized:

            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
        """
        from dendropy.utility.iosys import require_format_from_kwargs
        from dendropy.dataio import write_tree_list
        write_tree_list(format=require_format_from_kwargs(kwargs), tree_list=self, **kwargs)

    def reindex_subcomponent_taxa(self):
        """
        Synchronizes `TaxonSet` of member trees to `taxon_set` of self.
        """
        for t in self:
            self.reindex_tree_taxa(t)

    def reindex_tree_taxa(self, tree):
        """
        Synchronizes `tree` TaxonSet with self.
        """
        ti_mutable = self.taxon_set._is_mutable
        self.taxon_set._is_mutable = True
        tree.reindex_taxa(taxon_set=self.taxon_set, clear=False)
        self.taxon_set._is_mutable = ti_mutable

    def __setitem__(self, key, tree):
        """
        Homogeneity of taxon domain specified by `self.taxon_set`.
        """
        tree.reindex_tree_taxa(tree)
        list.__setitem__(self, key, tree)

    def extend(self, tree_list, reindex_taxa=True):
        """
        Homogeneity of taxon domain specified by `self.taxon_set`.
        """
        for t in tree_list:
            if reindex_taxa:
                self.reindex_tree_taxa(t)
            self.append(t)

    def append(self, tree, reindex_taxa=True):
        """
        Homogeneity of taxon domain specified by `self.taxon_set`.
        """
        if reindex_taxa:
            self.reindex_tree_taxa(tree)
        self[len(self):] = [tree]

    def as_python_source(self, tree_list_name=None, tree_list_args=None, oids=False):
        """
        Returns string that will rebuild this tree list in Python.
        """
        p = []

        if tree_list_name is None:
            tree_list_name = "tree_list_%s" % id(self)

        p.append("%s = dendropy.TreeList(label=%s%s%s)" \
            % (tree_list_name,
               ('"' + self.label +'"') if self.label is not None else "None",
               (', oid="%s"' % self.oid) if oids else "",
               (", " + tree_list_args) if tree_list_args is not None else ""))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        taxon_map = {}
        for taxon in self.taxon_set:
            tobj_name = taxon_obj_namer(taxon)
            p.append("%s = %s.taxon_set.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_list_name,
                   ('"' + taxon.label +'"') if taxon.label is not None else "None",
                   (', oid="%s"' % taxon.oid) if oids else ""))
            taxon_map[taxon] = tobj_name

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for tree in self:
            tree_obj_name = "tree_%s" % id(tree)
            p.append("%s = dendropy.Tree(label=%s, taxon_set=%s.taxon_set%s)" \
                % (tree_obj_name,
                   ('"' + tree.label +'"') if tree.label is not None else "None",
                   tree_list_name,
                   (', oid="%s"' % tree.oid) if oids else ""))
            p.append("%s.append(%s, reindex_taxa=False)" % (tree_list_name, tree_obj_name))
            if oids:
                p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, tree.seed_node.oid))
            for node in tree.preorder_node_iter():
                for child in node.child_nodes():
                    p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                            (node_obj_namer(child),
                             ("%s.seed_node" % tree_obj_name) if node is tree.seed_node else node_obj_namer(node),
                             ('"' + child.label +'"') if child.label is not None else "None",
                             taxon_obj_namer(child.taxon) if child.taxon is not None else "None",
                             child.edge.length,
                         (', oid="%s"' % child.oid) if oids else ""))
                    if oids:
                        p.append('%s.edge.oid = "%s"' % (node_obj_namer(child), child.edge.oid))

        return "\n".join(p)

##############################################################################
## Tree

class Tree(TaxonSetLinked, iosys.Readable, iosys.Writeable):
    """
    Fundamental class that encapsulates functionality and
    attributes need for working with trees.  A `Tree` contains a
    `seed_node` attribute (from which the entire tree springs), which
    may or may not be the root node. The distinction is not
    consequential in the current implementation, which identifies the
    root node as a node without `child_node` objects.
    """

    ###########################################################################
    ## Static methods

    def ancestor(node1, node2):
        """
        Returns the most-recent common ancestor node of node1 and
        node2.
        """
        mrca_node = None
        for node1_anc in Node.ancestor_iter(node1, inclusive=True):
            for node2_anc in Node.ancestor_iter(node2, inclusive=True):
                if node1_anc == node2_anc:
                    return node1_anc
        return None

    ancestor = staticmethod(ancestor)

    ###########################################################################
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        Initializes a new Tree object, optionally constructing it by cloning
        another Tree object if this is passed as the first argument, or
        out of a data source if `stream` and `format` are keyword arguments are
        passed with a file-like object and a format-specification string object
        values respectively.

        If `stream` and `format` keyword arguments are given, will
        construct this `Tree` object from `format`-formatted source
        given by file-like object `stream`. `format` must be a
        recognized and tree file format, such as `nexus`, `newick`, etc,
        for which a specialized tree list writer is available. If this
        is not implemented for the format specified, then a
        `UnsupportedFormatError` is raised. Other keywords will be
        passed to the underlying tree parser.

        Tree objects can thus be instantiated in the following ways::

            # /usr/bin/env python

            import StringIO
            import dendropy

            # empty tree
            t1 = Tree()

            # tree from data source
            t2 = Tree(stream=StringIO("((A,B),(C,D));"), format="newick")

            # passing keywords to underlying tree parser
            t3 = Tree(stream=StringIO("((A,B),(C,D));"),
                      format="newick",
                      taxon_set=t3.taxon_set,
                      encode_splits=True)

            # tree structure deep-copied from another tree
            t4 = Tree(t3)
            assert t4.taxon_set == t3.taxon_set # True: taxa are not deep-copied
            assert t4.oid != t3.oid # True: oid's will be different

            # can also call `read()` on a Tree object
            t5 = Tree()
            t5.read(StringIO("((A,B),(C,D));"), "newick")
            t5.read_from_string("((A,B),(C,D));", "newick")
            t5.read_from_path("mle.tre", "newick")

        """
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        iosys.Writeable.__init__(self)
        iosys.Readable.__init__(self)
        self.seed_node = Node(edge=Edge())
        self.length_type = None
        self.is_rooted = False

        if len(args) > 1:
            raise errors.TooManyArgumentsError(self.__class__.__name__, 1, args)
        if len(args) == 1:
            if "stream" in kwargs or "format" in kwargs:
                raise errors.MultipleInitializationSourceError(self.__class__.__name__, args[0])
            if isinstance(args[0], Node):
                self.seed_node = args[0]
            elif isinstance(args[0], Tree):
                self.clone_from(args[0])
            else:
                raise errors.InvalidArgumentTypeError(self.__class__.__name__, args[0])
        else:
            self.process_source_kwargs(**kwargs)
        if "oid" in kwargs:
            self.oid = kwargs["oid"]
        if "label" in kwargs:
            self.label = kwargs["label"]

    ###########################################################################
    ## I/O and Representation

    def clone_from(self, other):
        """
        Clones the structure and properties of `Tree` object `other`.
        """
        t = copy.deepcopy(other)
        self.__dict__ = t.__dict__
        return self

    def __deepcopy__(self, memo):
        # we treat the taxa as immutable and copy the reference even in a deepcopy
        o = self.__class__(taxon_set=self.taxon_set)
        memo[id(self)] = o
        memo[id(self.taxon_set)] = o.taxon_set
        for i, t in enumerate(self.taxon_set):
            memo[id(t)] = o.taxon_set[i]
#        if self.seed_node is not None:
#            new_v = copy.deepcopy(self.seed_node, memo)
#            o.seed_node = new_v
#        else:
#            o.seed_node = None
#        memo[id(self.seed_node)] = o.seed_node
        for k, v in self.__dict__.iteritems():
            if k not in ['taxon_set', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

    def __str__(self):
        "Dump Newick string."
        return self.as_newick_str()

    def read(self, stream, format, **kwargs):
        """
        Populates/constructs objects of this type from `format`-formatted
        data in the file-like object source `stream`.

        Recognized keywords arguments are:

            - `taxon_set` specifies the `TaxonSet` object to be attached to the
               trees parsed and manage their taxa. If not specified, then the
               `TaxonSet` object currently associated with the tree will be used.
            - `encode_splits` specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `translate_dict` should provide a dictionary mapping taxon numbers (as
               found in the source) to taxon labels (as defined in the source).
            - `rooted` specifies the default rooting interpretation of the tree (see
               `dendropy.dataio.nexustokenizer` for details).
            - `finish_node_func` is a function that will be applied to each node
               after it has been constructed.
            - `edge_len_type` specifies the type of the edge lengths (int or float)

        If the source defines multiple trees, only the first one will be
        returned unless the keyword `index` is used to specify the
        0-based index of the tree to be returned. If `index` >= number
        of trees, a KeyError is raised.
        """
        from dendropy.utility import iosys
        from dendropy.dataio import tree_source_iter
        if "index" in kwargs:
            index = kwargs.get("index")
            del(kwargs["index"])
        else:
            index = 0
        if "taxon_set" not in kwargs:
            kwargs["taxon_set"] = self.taxon_set
        else:
            self.taxon_set = kwargs["taxon_set"]
        titer = tree_source_iter(stream=stream, format=format, **kwargs)
        count = 0
        t = None
        while count <= index:
            try:
                t = titer.next()
            except StopIteration:
                raise KeyError("0-based index out of bounds: %d (trees=%d, index=[0, %d])" % (index, count, count-1))
            else:
                count += 1
        self.__dict__ = t.__dict__
        return self

    def write(self, **kwargs):
        """
        Writes out `Tree` in `format` to a destination described by
        one of: `file` or `path`:

            - `file`: A file- or file-like object.
            - `path`: A string specifying the path to a file.

        `format` must be a recognized and tree file format, such as `nexus`,
        `newick`, etc, for which a specialized tree list writer is
        available. If this is not implemented for the format specified, then
        a `UnsupportedFormatError` is raised.

        Additionally, for some formats, the following keywords are recognized:

            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
        """
        from dendropy.utility.iosys import require_format_from_kwargs
        from dendropy.dataio import write_tree_list
        tree_list = TreeList(taxon_set=self.taxon_set)
        tree_list.append(self)
        write_tree_list(format=require_format_from_kwargs(kwargs), tree_list=tree_list, **kwargs)

    ###########################################################################
    ## Getting/accessing methods

    def nodes(self, cmp_fn=None, filter_fn=None):
        "Returns list of nodes on the tree, sorted using cmp_fn."
        nodes = [node for node in self.preorder_node_iter(filter_fn)]
        if cmp_fn:
            nodes.sort(cmp_fn)
        return nodes

    def leaf_nodes(self):
        "Returns list of leaf_nodes on the tree."
        return [leaf for leaf in self.leaf_iter()]

    def internal_nodes(self):
        "Returns list of internal node in the tree."
        return self.nodes(filter_fn=lambda x : not x.is_leaf())

    def find_node_for_taxon(self, taxon):
        for node in self.preorder_node_iter():
            try:
                if node.taxon is taxon:
                    return node
            except:
                pass
        return None

    def find_node(self, filter_fn):
        """
        Finds the first node for which filter_fn(node) = True.
        For example, if::

            filter_fn = lambda n: hasattr(n, 'genes') and n.genes is not None

        then::

            t.find_node(filter_fn=filter_fn)

        will return all nodes which have an attributed 'genes' and this value
        is not None.
        """
        for node in self.preorder_node_iter(filter_fn):
            return node
        return None

    def find_taxon_node(self, taxon_filter_fn=None):
        "Finds the first node for which taxon_filter_fn(node.taxon) == True."
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None

    def find_edge(self, oid):
        "Finds the first edge with matching id."
        for edge in self.preorder_edge_iter():
            if edge.oid == oid:
                return edge
        return None

    def get_edge_set(self, filter_fn=None):
        """Returns the set of edges that are currently in the tree.

        Note: the returned set acts like a shallow copy of the edge set (adding
        or deleting elements from the set does not change the tree, but
        modifying the elements does).
        """
        return set([i in self.preorder_edge_iter(filter_fn=filter_fn)])

    def get_node_set(self, filter_fn=None):
        """Returns the set of nodes that are currently in the tree

        Note: the returned set acts like a shallow copy of the edge set (adding
        or deleting elements from the set does not change the tree, but
        modifying the elements does).
        """
        return set([i in self.preorder_node_iter(filter_fn=filter_fn)])

    ###########################################################################
    ## Node iterators

    def __iter__(self):
        return self.postorder_node_iter()

    def preorder_node_iter(self, filter_fn=None):
        "Returns preorder iterator over tree nodes."
        for node in self.seed_node.preorder_iter(self.seed_node, filter_fn):
            yield node

    def postorder_node_iter(self, filter_fn=None):
        "Returns postorder iterator over tree nodes."
        for node in self.seed_node.postorder_iter(self.seed_node, filter_fn):
            yield node

    def level_order_node_iter(self, filter_fn=None):
        "Returns level-order iterator over tree nodes."
        for node in self.seed_node.level_order_iter(self.seed_node, filter_fn):
            yield node

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over tree leaf_nodes (order determined by
        postorder tree-traversal).
        """
        for node in self.seed_node.leaf_iter(self.seed_node, filter_fn):
            yield node

    ###########################################################################
    ## Edge iterators

    def preorder_edge_iter(self, filter_fn=None):
        "Returns preorder iterator over tree edges."
        for node in self.seed_node.preorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def postorder_edge_iter(self, filter_fn=None):
        "Returns postorder iterator over tree edges."
        for node in self.seed_node.postorder_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def level_order_edge_iter(self, filter_fn=None):
        "Returns level-order iterator over tree edges."
        for node in self.seed_node.level_order_iter(self.seed_node):
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    ###########################################################################
    ## Information/Utilities

    def add_ages_to_nodes(self, attr_name='age', check_prec=0.0000001):
        """
        Takes an ultrametric `tree` and adds a attribute named `attr` to
        each node, with the value equal to the sum of edge lengths from the
        node to the tips. If the lengths of different paths to the node
        differ by more than `check_prec`, then a ValueError exception
        will be raised indicating deviation from ultrametricity. If
        `check_prec` is negative or False, then this check will be
        skipped.
        """
        for node in self.postorder_node_iter():
            ch = node.child_nodes()
            if len(ch) == 0:
                setattr(node, attr_name, 0.0)
            else:
                first_child = ch[0]
                setattr(node, attr_name, getattr(first_child, attr_name) + first_child.edge.length)
                if not (check_prec < 0 or check_prec == False):
                    for nnd in ch[1:]:
                        ocnd = getattr(nnd, attr_name) + nnd.edge.length
                        if abs(getattr(node, attr_name) - ocnd) > check_prec:
                            raise ValueError("Tree is not ultrametric")

    ###########################################################################
    ## Taxa Management

    def infer_taxa(self):
        """
        Returns a new TaxonSet object populated with taxa from this
        tree.
        """
        taxon_set = taxon.TaxonSet()
        for node in self.postorder_node_iter():
            if node.taxon is not None and (node.taxon not in taxon_set):
                taxon_set.add(node.taxon)
        self.taxon_set = taxon_set
        return taxon_set

    def reindex_subcomponent_taxa(self):
        """
        Reassigns node taxon objects
        """
        for node in self.postorder_node_iter():
            t = node.taxon
            if t:
                node.taxon = self.taxon_set.require_taxon(label=t.label)

    ###########################################################################
    ## Structure

    def deroot(self):
        "Converts a degree-2 node at the root to a degree-3 node."
        seed_node = self.seed_node
        if not seed_node:
            return
        child_nodes = seed_node.child_nodes()
        if len(child_nodes) != 2:
            return

        if len(child_nodes[1].child_nodes()) >= 2:
            to_keep, to_del = child_nodes
        elif len(child_nodes[0].child_nodes()) >= 2:
            to_del, to_keep = child_nodes
        else:
            return
        to_del_edge = to_del.edge
        try:
            to_keep.edge.length += to_del_edge.length
        except:
            pass
        from dendropy.treemanip import collapse_edge
        collapse_edge(to_del_edge)

    def reroot_at(self, nd, splits=False, delete_deg_two=True):
        """Takes an internal node, `nd` that must already be in the tree and
        reroots the tree such that `nd` is the `seed_node` of the tree.

        If `splits` is True, then the edges' `clade_mask` and the tree's
            `split_edges` attributes will be updated."""
        old_par = nd.parent_node
        if old_par is None:
            return
        if splits:
            taxa_mask = self.seed_node.edge.clade_mask
        to_edge_dict = None
        if splits:
            to_edge_dict = getattr(self, "split_edges", None)

        if old_par is self.seed_node:
            root_children = old_par.child_nodes()
            if len(root_children) == 2 and delete_deg_two:
                # root (old_par) was of degree 2, thus we need to suppress the
                #   node
                fc = root_children[0]
                if fc is nd:
                    sister = root_children[1]
                else:
                    assert root_children[1] is nd
                    sister = fc
                if nd.edge.length:
                    sister.edge.length += nd.edge.length
                edge_to_del = nd.edge
                nd.edge = old_par.edge
                if splits:
                    assert nd.edge.clade_mask == taxa_mask
                if to_edge_dict:
                    del to_edge_dict[edge_to_del.clade_mask]
                nd.add_child(sister, edge_length=sister.edge.length)
                self.seed_node = nd
                return
        else:
            self.reroot_at(old_par, splits=splits, delete_deg_two=delete_deg_two)
        old_par.edge, nd.edge = nd.edge, old_par.edge
        e = old_par.edge
        if splits:
            if to_edge_dict:
                del to_edge_dict[e.clade_mask]
            e.clade_mask = (~(e.clade_mask)) & taxa_mask
            if to_edge_dict:
                to_edge_dict[e.clade_mask] = e
            assert nd.edge.clade_mask == taxa_mask
        old_par.remove_child(nd)
        nd.add_child(old_par, edge_length=e.length)
        self.seed_node = nd

    def to_outgroup_position(self, nd, splits=False, delete_deg_two=True):
        """Reroots the tree at the parend of `nd` and makes `nd` the first child
        of the new root.  This is just a convenience function to make it easy
        to place a clade as the first child under the root.

        Assumes that `nd` and `nd.parent_node` and are in the tree

        If `splits` is True, then the edges' `clade_mask` and the tree's
            `split_edges` attributes will be updated.
        If `delete_deg_two` is True and the old root of the tree has an
            outdegree of 2, then the node will be removed from the tree.
        """
        p = nd.parent_node
        assert p is not None
        self.reroot_at(p, splits=splits)
        p.remove_child(nd)
        p.add_child(nd, edge_length=nd.edge.length, pos=0)

    ###########################################################################
    ## For debugging

    def as_newick_str(self, **kwargs):
        """kwargs["reverse_translate"] can be function that takes a taxon and
           returns the label to appear in the tree."""
        return self.seed_node.as_newick_str(**kwargs)

    def assign_node_labels_from_taxon_or_oid(self):
        for nd in self.postorder_node_iter():
            if nd.label is not None:
                continue
            if nd.taxon is not None:
                nd.label = nd.taxon.label
            else:
                nd.label = nd.oid

    def as_python_source(self, tree_obj_name=None, tree_args=None, oids=False):
        """
        Returns string that will rebuild this tree in Python.
        """
        p = []

        if tree_obj_name is None:
            tree_obj_name = "tree_%s" % id(self)
        p.append("%s = dendropy.Tree(label=%s%s)" \
            % (tree_obj_name,
               ('"' + self.label +'"') if self.label is not None else "None",
               (', oid="%s"' % self.oid) if oids else "",
               (", " + tree_args) if tree_args is not None else ""))
        if oids:
            p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, self.seed_node.oid))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        for taxon in self.taxon_set:
            tobj_name = taxon_obj_namer(taxon)
            p.append("%s = %s.taxon_set.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_obj_name,
                   ('"' + taxon.label +'"') if taxon.label is not None else "None",
                   (', oid="%s"' % taxon.oid) if oids else ""))

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for node in self.preorder_node_iter():
            for child in node.child_nodes():
                p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                        (node_obj_namer(child),
                         ("%s.seed_node" % tree_obj_name) if node is self.seed_node else node_obj_namer(node),
                         ('"' + node.label +'"') if child.label is not None else "None",
                         taxon_obj_namer(child.taxon) if child.taxon is not None else "None",
                         child.edge.length,
                         (', oid="%s"' % child.oid) if oids else ""))
                if oids:
                    p.append('%s.edge.oid = "%s"' % (node_obj_namer(child), child.edge.oid))

        return "\n".join(p)

    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        if kwargs.get("splits"):
            if not kwargs.get("taxon_set"):
                kwargs["taxon_set"] = self.taxon_set
        self.seed_node.write_indented_form(out, **kwargs)

    def debug_check_tree(self, logger_obj=None, **kwargs):
        import logging, inspect
        if logger_obj and logger_obj.isEnabledFor(logging.DEBUG):
            try:
                assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)
            except:
                calling_frame = inspect.currentframe().f_back
                co = calling_frame.f_code
                emsg = "\nCalled from file %s, line %d, in %s" % (co.co_filename, calling_frame.f_lineno, co.co_name)
                _LOG.debug("%s" % str(self))
                _LOG.debug("%s" % self.get_indented_form(**kwargs))
        assert self._debug_tree_is_valid(logger_obj=logger_obj, **kwargs)

    def _debug_tree_is_valid(self, **kwargs):
        """Performs sanity-checks of the tree data structure.

        kwargs:
            `splits` if True specifies that the split_edge and clade_mask attributes
                are checked.
        """
        check_splits = kwargs.get('splits', False)
        taxon_set = kwargs.get('taxon_set')
        if taxon_set is None:
            taxon_set = self.taxon_set
        if check_splits:
            taxa_mask = self.seed_node.edge.clade_mask
        nodes = set()
        edges = set()
        curr_node = self.seed_node
        assert(curr_node.parent_node is None)
        assert(curr_node.edge.tail_node is None)
        ancestors = []
        siblings = []
        while curr_node:
            curr_edge = curr_node.edge
            assert(curr_edge not in edges)
            edges.add(curr_edge)
            assert(curr_node not in nodes)
            nodes.add(curr_node)
            assert(curr_edge.tail_node is curr_node.parent_node)
            assert(curr_edge.head_node is curr_node)
            if check_splits:
                cm = 0
                clade_mask = curr_edge.clade_mask
                assert((clade_mask | taxa_mask) == taxa_mask)
            c = curr_node.child_nodes()
            if c:
                for child in c:
                    assert child.parent_node is curr_node
                    if check_splits:
                        cm |= child.edge.clade_mask
            elif check_splits:
                cm = taxon_set.taxon_bitmask(curr_node.taxon)
            if check_splits:
                assert((cm & taxa_mask) == clade_mask)
                assert self.split_edges[clade_mask] == curr_edge
            curr_node, level = _preorder_list_manip(curr_node, siblings, ancestors)
        if check_splits:
            for s, e in self.split_edges.iteritems():
                assert(e in edges)
        return True

##############################################################################
## Node

class Node(TaxonLinked):
    """
    A node on a tree, implementing only fundamental behaviour and
    properties.
    """

    ### ITERATORS ############################################################

    def preorder_iter(node, filter_fn=None):
        """
        Preorder traversal of the node and its child_nodes.  Returns node
        and all descendants such that node is returned before node's
        child_nodes (and their child_nodes). Filtered by filter_fn: node is
        only returned if no filter_fn is given or if filter_fn returns
        True.
        """
        if not node:
            return
        stack = [node]
        while stack:
            node = stack.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            child_nodes.extend(stack)
            stack = child_nodes
    preorder_iter = staticmethod(preorder_iter)

    def postorder_iter(node, filter_fn=None):
        """
        Postorder traversal of the node and its child_nodes.  Returns node
        and all descendants such that node's child_nodes (and their
        child_nodes) are visited before node.  Filtered by filter_fn:
        node is only returned if no filter_fn is given or if filter_fn
        returns True.
        """
        stack = [(node, False)]
        while stack:
            node, state = stack.pop(0)
            if state:
                if filter_fn is None or filter_fn(node):
                    yield node
            else:
                stack.insert(0, (node, True))
                child_nodes = [(n, False) for n in node.child_nodes()]
                child_nodes.extend(stack)
                stack = child_nodes
    postorder_iter = staticmethod(postorder_iter)

    def leaf_iter(start_nd, filter_fn=None):
        """
        Returns an iterator over the leaf_nodes that are descendants `of start_nd`
        (order determined by postorder tree-traversal).
        """
        if filter_fn:
            filter_fn = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            filter_fn = lambda x: x.is_leaf() and x or None
        for node in start_nd.postorder_iter(start_nd, filter_fn):
            yield node

    leaf_iter = staticmethod(leaf_iter)

    def level_order_iter(node, filter_fn=None):
        """
        Level-order traversal of the node and its child_nodes. Filtered
        by filter_fn: node is only returned if no filter_fn is given
        or if filter_fn returns True
        """
        if filter_fn is None or filter_fn(node):
            yield node
        remaining = node.child_nodes()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            remaining.extend(child_nodes)

    level_order_iter = staticmethod(level_order_iter)

    def ancestor_iter(node, filter_fn=None, inclusive=True):
        """
        Returns all ancestors of node. If `inclusive` is True, `node`
        is returned as the first item of the sequence.
        """
        if inclusive:
            yield node
        while node is not None:
            node = node.parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    ancestor_iter = staticmethod(ancestor_iter)

    ## UTILITIES #############################################################

    def nodeset_hash(nodes, attribute='oid'):
        """
        Returns a hash of a set of nodes, based on the given
        attribute.
        """
        tags = []
        for node in nodes:
            if hasattr(node, attribute) and getattr(node, attribute) != None:
                value = getattr(node, attribute)
                tags.append(str(value))
        tags.sort()
        return '+'.join(tags)

    nodeset_hash = staticmethod(nodeset_hash)

    ## INSTANCE METHODS########################################################

    def __init__(self, **kwargs):
        TaxonLinked.__init__(self,
                             taxon=kwargs.get("taxon", None),
                             label=kwargs.get("label", None),
                             oid=kwargs.get("oid", None))
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        self.edge = kwargs.get("edge", Edge(head_node=self))
        self._edge.head_node = self

    def __deepcopy__(self, memo):
        o = self.__class__(taxon=self.taxon)
        memo[id(self)] = o
        if self.taxon is not None:
            memo[id(self.taxon)] = o.taxon
        for k, v in self.__dict__.iteritems():
            if not k in ['_child_nodes', '_taxon', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        for c in self.child_nodes():
            o.add_child(copy.deepcopy(c, memo))
        memo[id(self._child_nodes)] = o._child_nodes
        memo[id(self._oid)] = o._oid
        return o

    def collapse_neighborhood(self, dist):
        if dist < 1:
            return
        children = self.child_nodes()
        for ch in children:
            if not ch.is_leaf():
                ch.edge.collapse()
        if self.parent_node:
            p = self.parent_node
            self.edge.collapse()
            p.collapse_neighborhood(dist -1)
        else:
            self.collapse_neighborhood(dist - 1)

    def is_leaf(self):
        "Returns True if the node has no child_nodes"
        return bool(not self._child_nodes)

    def is_internal(self):
        "Returns True if the node has child_nodes"
        return bool(self._child_nodes)

    ## Low-level methods for manipulating structure ##

    def _get_edge(self):
        "Returns the edge subtending this node."
        return self._edge

    def _set_edge(self, edge=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self._edge = edge
        if edge:
            edge.head_node = self

    def _get_edge_length(self):
        "Returns the length of the edge  subtending this node."
        return self._edge.length

    def _set_edge_length(self, v=None):
        """
        Sets the edge subtending this node, and sets head_node of
        `edge` to point to self.
        """
        self._edge.length = v

    edge = property(_get_edge, _set_edge)
    edge_length = property(_get_edge_length, _set_edge_length)

    def child_nodes(self):
        "Returns the a shallow-copy list of all child nodes."
        return list(self._child_nodes)

    def set_children(self, child_nodes):
        """
        Sets the child_nodes for this node.
        Side effects:
            - sets the parent of each child node to this node
            - sets the tail node of each child to self
        """
        self._child_nodes = child_nodes
        for nidx in range(len(self._child_nodes)):
            self._child_nodes[nidx].parent_node = self
            self._child_nodes[nidx].edge.tail_node = self

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self._parent_node

    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self._parent_node = parent
        self.edge.tail_node = parent

    parent_node = property(_get_parent_node, _set_parent_node)

    def get_incident_edges(self):
        e = [c.edge for c in self._child_nodes]
        e.append(self.edge)
        return e

    def get_adjacent_nodes(self):
        n = [c for c in self._child_nodes]
        if self.parent_node:
            n.append(self.parent_node)
        return n

    def add_child(self, node, edge_length=None, pos=None):
        """
        Adds a child node to this node. Results in the parent_node and
        containing_tree of the node being attached set to this node.
        If `edge_length` is given, then the new child's edge length is
        set to this. Returns node that was just attached.
        """
        node.parent_node = self
        if edge_length != None:
            node.edge_length = edge_length
        if pos is None:
            self._child_nodes.append(node)
        else:
            self._child_nodes.insert(pos, node)
        return node

    def new_child(self, **kwargs):
        """
        Convenience class to create and add a new child to this node. Keyword
        arguments `label`, `oid` and `taxon` will be passed to Node().
        `edge_length`, if given, will specify the length of the subtending
        edge of this child.
        """
        node = self.__class__(**kwargs)
        return self.add_child(node, edge_length=kwargs.get("edge_length", None))

    def remove_child(self, node, suppress_deg_two=False):
        """
        Removes a node from this nodes child set. Results in the
        parent of the node being removed set to None.

        Returns the node removed.

        `suppress_deg_two` should only be called on unrooted trees.
        """
        if not node:
            raise Exception("Tried to remove an non-existing or null node")
        children = self._child_nodes
        if node in children:
            node.parent_node = None
            node.edge.tail_node = None
            index = children.index(node)
#             if index > 0:
#                 self._child_nodes[index-1].next_sib = None
            children.remove(node)
            if suppress_deg_two:
                if self.parent_node:
                    if len(children) == 1:
                        child = children[0]
                        pos = self.parent_node._child_nodes.index(self)
                        self.parent_node.add_child(child, pos=pos)
                        self.parent_node.remove_child(self, suppress_deg_two=False)
                        try:
                            child.edge.length += self.edge.length
                        except:
                            pass
                        self._child_nodes = []
                else:
                    to_remove = None
                    if len(children) == 2:
                        if children[0].is_internal():
                            to_remove = children[0]
                            other = children[1]
                        elif children[1].is_internal():
                            to_remove = children[1]
                            other = children[0]
                    if to_remove is not None:
                        try:
                            other.edge.length += to_remove.edge.length
                        except:
                            pass
                        pos = self._child_nodes.index(to_remove)
                        self.remove_child(to_remove, suppress_deg_two=False)
                        tr_children = to_remove._child_nodes
                        tr_children.reverse()
                        for c in tr_children:
                            self.add_child(c, pos=pos)
                        to_remove._child_nodes = []
        else:
            raise Exception("Tried to remove a node that is not listed as a child")
        return node

    def reversible_remove_child(self, node, suppress_deg_two=False):
        """
        This function is a (less-efficient) version of remove_child that also
        returns the data needed by reinsert_nodes to "undo" the removal.

        Returns a list of tuples.  The first element of each tuple is the
        node removed, the other elements are the information needed by
        `reinsert_nodes' in order to restore the tree to the same topology as
        it was before the call to `remove_child.` If `suppress_deg_two` is False
        then the returned list will contain only one item.

        `suppress_deg_two` should only be called on unrooted trees.
        """
        if not node:
            raise Exception("Tried to remove an non-existing or null node")
        children = self._child_nodes
        try:
            pos = children.index(node)
        except:
            raise Exception("Tried to remove a node that is not listed as a child")

        removed = [(node, self, pos, [], None)]
        node.parent_node = None
        node.edge.tail_node = None
#             if index > 0:
#                 self._child_nodes[index-1].next_sib = None
        children.remove(node)
        if suppress_deg_two:
            p = self.parent_node
            if p:
                if len(children) == 1:
                    child = children[0]
                    pos = p._child_nodes.index(self)
                    p.add_child(child, pos=pos)
                    self._child_nodes = []
                    p.remove_child(self, suppress_deg_two=False)
                    e = child.edge
                    try:
                        e.length += self.edge.length
                    except:
                        e = None
                    t = (self, p, pos, [child], e)
                    removed.append(t)
            else:
                to_remove = None
                if len(children) == 2:
                    if children[0].is_internal():
                        to_remove = children[0]
                        other = children[1]
                    elif children[1].is_internal():
                        to_remove = children[1]
                        other = children[0]
                if to_remove is not None:
                    e = other.edge
                    try:
                        e.length += to_remove.edge.length
                    except:
                        e = None
                    pos = self._child_nodes.index(to_remove)
                    self.remove_child(to_remove, suppress_deg_two=False)
                    tr_children = to_remove._child_nodes
                    to_remove._child_nodes = []
                    for n, c in enumerate(tr_children):
                        new_pos = pos + n
                        self.add_child(c, pos=new_pos)
                    t = (to_remove, self, pos, tr_children, e)
                    removed.append(t)


        return removed

    def reinsert_nodes(self, nd_connection_list):
        """T
        his function should be used to "undo" the effects of
        Node.reversible_remove_child NOTE: the behavior is only
        guaranteed if the tree has not been modified between the
        remove_child and reinsert_nodes calls! (or the tree has been
        restored such that the node/edge identities are identical to the
        state before the remove_child call.

        The order of info in each tuple is:

            0 - node removed
            1 - parent of node removed
            2 - pos in parent array
            3 - children of node removed that were "stolen"
            4 - edge that was lengthened by "stealing" length from node's edge
        """
        # we unroll the stack of operations
        for blob in nd_connection_list[-1::-1]:
            #_LOG.debug(blob)
            n, p, pos, children, e = blob
            for c in children:
                cp = c.parent_node
                if cp:
                    cp.remove_child(c)
                n.add_child(c)
            p.add_child(n, pos=pos)
            if e is not None:
                e.length -= n.edge.length

    ## Basic node metrics ##

    def distance_from_tip(self):
        """
        Sum of edge lengths from tip to node. If tree is not ultrametric
        (i.e., descendent edges have different lengths), then count the
        maximum of edge lengths. Note that the 'add_ages_to_nodes()' method
        of dendropy.trees.Tree() is a more efficient way of doing this over
        the whole tree.
        """
        if not self._child_nodes:
            return 0.0
        else:
            distance_from_tips = []
            for ch in self._child_nodes:
                if ch.edge.length is not None:
                    curr_edge_length = ch.edge_length
                else:
                    curr_edge_length = 0.0
                if not hasattr(ch, "_distance_from_tip"):
                    ch._distance_from_tip = ch.distance_from_tip()
                distance_from_tips.append(ch._distance_from_tip + curr_edge_length)
            self._distance_from_tip = float(max(distance_from_tips))
            return self._distance_from_tip

    def distance_from_root(self):
        """
        Sum of edge lengths from root. Right now, 'root' is taken to
        be a node with no parent node.
        """
        if self.parent_node and self.edge.length != None:
            if self.parent_node.distance_from_root == None:
                return float(self.edge.length)
            else:
                distance_from_root = float(self.edge.length)
                parent_node = self.parent_node
                # The root is identified when a node with no
                # parent is encountered. If we want to use some
                # other criteria (e.g., where a is_root property
                # is True), we modify it here.
                while parent_node:
                    if parent_node.edge.length != None:
                        distance_from_root = distance_from_root + float(parent_node.edge.length)
                    parent_node = parent_node.parent_node
                return distance_from_root
        elif not self.parent_node and self.edge.length != None:
            return float(self.edge.length)
        elif self.parent_node and self.edge.length == None:
            # what do we do here: parent node exists, but my
            # length does not?
            return float(self.parent_node.edge.length)
        elif not self.parent_node and self.edge.length == None:
            # no parent node, and no edge length
            return 0.0
        else:
            # WTF????
            return 0.0

    def level(self):
        "Number of nodes between self and root."
        if self.parent_node:
            return self.parent_node.level + 1
        else:
            return 0

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with self as the only member if self is a leaf).
        """
        return [node for node in \
                self.postorder_iter(self, \
                                    lambda x: bool(len(node.child_nodes)==0))]

    ###########################################################################
    ## For debugging we build-in a full-fledged NEWICK composition independent
    ## of the nexus/newick family of modules. Client code should prefer to
    ## use Newick/Nexus readers/writers, or Tree.write(), TreeList.write(),
    ## DataSet.write() etc.

    def as_newick_str(self, **kwargs):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        out = StringIO()
        self.write_newick(out, **kwargs)
        return out.getvalue()

    def write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK
        statement according to the given formatting rules.
        """
        edge_lengths = kwargs.get('edge_lengths', True)
        child_nodes = self.child_nodes()
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                child.write_newick(out, **kwargs)
            out.write(')')

        out.write(self.get_node_str(**kwargs))
        if edge_lengths:
            e = self.edge
            if e:
                sel = e.length
                if sel is not None:
                    fmt = kwargs.get('edge_length_formatter', None)
                    if fmt:
                        out.write(":%s" % fmt(sel))
                    else:
                        s = ""
                        try:
                            s = float(sel)
                            s = str(s)
                        except ValueError:
                            s = str(sel)
                        if s:
                            out.write(":%s" % s)

    def get_node_str(self, **kwargs):
        """returns a string that is an identifier for the node.  This is called
        by the newick-writing functions, so the kwargs that affect how node
        labels show up in a newick string are the same ones used here:
        `include_internal_labels` is a Boolean.
        """
        is_leaf = (len(self._child_nodes) == 0)
        include_internal_labels = kwargs.get("include_internal_labels")
        if (not is_leaf) and (not include_internal_labels):
            return ""
        try:
            t = self.taxon
            rt = kwargs.get("reverse_translate")
            if rt:
                tag = rt(t)
            else:
                tag = t.label

        except AttributeError:
            tag = ""
            try:
                tag = self.label
            except AttributeError:
                if not is_leaf:
                    tag = self.oid
        if "raw_labels" in kwargs:
            return tag
        return texttools.escape_nexus_token(tag)

    ###########################################################################
    ## alternate representation of tree structure for debugging

    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        clade_masks = kwargs.get("splits", True)
        level = kwargs.get("level", 0)
        ancestors = []
        siblings = []
        n = self
        while n is not None:
            n._write_indented_form_line(out, level, **kwargs)
            n, lev = _preorder_list_manip(n, siblings, ancestors)
            level += lev

    def _get_indented_form_line(self, level, **kwargs):
        out = StringIO()
        self._write_indented_form_line(out, level, **kwargs)
        return out.getvalue()

    def _write_indented_form_line(self, out, level, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        label = format_node(self, **kwargs)
        if kwargs.get("splits"):
            cm = "%s " % format_split(self.edge.clade_mask, **kwargs)
        else:
            cm = ""
        out.write("%s%s%s\n" % ( cm, indentation*level, label))
##############################################################################
## Edge

class Edge(IdTagged):
    """
    An edge on a tree. This class implements only the core
    functionality needed for trees.
    """

    ## CLASS METHODS  ########################################################

    def __init__(self, **kwargs):
        """
        Creates an edge from tail_node to head_node.  Modified from
        arbol.
        """
        IdTagged.__init__(self, oid=kwargs.get("oid", None), label=kwargs.get("label", None))
        self.tail_node = kwargs.get("tail_node", None)
        self.head_node = kwargs.get("head_node", None)
        self.rootedge = kwargs.get("rootedge", False)
        self.length = kwargs.get("length", None)

    def __deepcopy__(self, memo):
        o = self.__class__()
        memo[id(self)] = o
        o.tail_node = copy.deepcopy(self.tail_node, memo)
        o.head_node = copy.deepcopy(self.head_node, memo)
        o.length = copy.deepcopy(self.length, memo)
        o.rootedge = copy.deepcopy(self.rootedge, memo)
        memo[id(self.tail_node)] = o.tail_node
        memo[id(self.head_node)] = o.head_node
        memo[id(self.length)] = o.length
        memo[id(self.rootedge)] = o.rootedge
        for k, v in self.__dict__.iteritems():
            if not k in ['tail_node', 'head_node', 'length', 'rootedge', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

    def collapse(self):
        h = self.head_node
        if h.is_leaf():
            return
        t = self.tail_node
        if t is None:
            return
        c = h.child_nodes()
        pc = t.child_nodes()
        pos = len(pc)
        try:
            pos = pc.index(h)
        except:
            pass
        for i, ch in enumerate(c):
            t.add_child(ch, pos=pos + i)
        t.remove_child(h)

    def new_edge(self, oid=None):
        "Returns a new edge object of the same class of this edge."
        edge = self.__class__()
        edge.oid = oid
        return edge

    def invert(self):
        self.head_node, self.tail_node = self.tail_node, self.head_node

    def is_terminal(self):
        "Returns True if the head node has no children"
        return bool(self.head_node and self.head_node.is_leaf())

    def is_internal(self):
        "Returns True if the head node has children"
        return bool(self.head_node and not self.head_node.is_leaf())

    def get_adjacent_edges(self):
        'Returns a list of all edges that "share" a node with `self`'
        he = [i for i in self.head_node.get_incident_edges() if i is not self]
        te = [i for i in self.tail_node.get_incident_edges() if i is not self]
        he.extend(te)
        return he
    adjacent_edges = property(get_adjacent_edges)

def _preorder_list_manip(n, siblings, ancestors):
    """
    Helper function for recursion free preorder traversal, that does
    not rely on attributes of the node other than child_nodes() (thus it
    is useful for debuggging).

    Returns the next node (or None) and the number of levels toward the
    root the function "moved".
    """
    levels_moved = 0
    c = n.child_nodes()
    if c:
        levels_moved += 1
        ancestors.append(list(siblings))
        del siblings[:]
        siblings.extend(c[1:])
        return c[0], levels_moved
    while not siblings:
        if ancestors:
            levels_moved -= 1
            del siblings[:]
            siblings.extend(ancestors.pop())
        else:
            return None, levels_moved
    return siblings.pop(0), levels_moved

def format_node(nd, **kwargs):
    if nd.is_leaf():
        t = nd.taxon
        if t:
            label = t.label
        else:
            label = "anonymous leaf"
    else:
        label = "* %s" % str(nd.oid)
    return label

def format_split(split, width=None, **kwargs):
    from dendropy.splitcalc import split_as_string
    if width is None:
        width = len(kwargs.get("taxon_set"))
    s = split_as_string(split, width, symbol1=kwargs.get("off_symbol"), symbol2=kwargs.get("on_symbol"))
    return s

