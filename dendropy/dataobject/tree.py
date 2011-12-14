#! /usr/bin/env pytho:

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
This module handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from cStringIO import StringIO
import copy
import bisect
import re

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

from dendropy.utility import GLOBAL_RNG
from dendropy.utility import iosys
from dendropy.utility import error
from dendropy.utility import textutils
from dendropy.utility import termutils
from dendropy.dataobject.base import IdTagged
from dendropy.dataobject.taxon import TaxonSet, TaxonSetLinked, TaxonLinked
from dendropy import treesplit

##############################################################################
## TreeList

class TreeList(list, TaxonSetLinked, iosys.Readable, iosys.Writeable):
    """
    Collects and coordinates a list of trees with the associated with the
    same set of taxa.
    """

    def __init__(self, *args, **kwargs):
        """
        __init__ creates a new TreeList object, populating it with any iterable
        container with Tree object members passed as unnamed argument, or
        from a data source if `stream` and `schema` are passed.

        If passed an iterable container, the objects in that container must be
        of type `Tree` (or derived). If the container is of type `TreeList`,
        then, because each `Tree` object must have the same `TaxonSet`
        reference as the containing `TreeList`, the trees in the container
        passed as an initialization argument will be **deep**-copied (except
        for associated TaxonSet and Taxon objects, which will be
        shallow-copied). If the container is any other type of iterable, then
        the `Tree` objects will be **shallow**-copied.

        TreeList objects can directly thus be instantiated in the
        following ways::

            # /usr/bin/env python

            import StringIO
            from dendropy import TaxonSet, Tree, TreeList

            # empty tree
            tlst1 = TreeList()

            # populated from list of Tree objects
            t1 = Tree(stream=StringIO("((A,B),(C,D))"), schema="newick")
            t2 = Tree(stream=StringIO("((A,C),(B,D))"), schema="newick")
            tlist2 = TreeList([t1, t2])

            # tree from data source
            tlst3 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"), schema="newick") # same

            # passing keywords to underlying tree parser
            tlst4 = TreeList(stream=StringIO("((A,B),(C,D));((A,C),(B,D));"),
                             schema="newick",
                             taxon_set=tlst3.taxon_set,
                             encode_splits=True)

            # deep-copied (but shallow-copy taxa) from another tree list
            tlst5 = TreeList(t4)

            # same
            tls6 = TreeList([Tree(t) for t in tlst5])

            # the canonical way to instantiate a TreeList from a data source
            # is the use the `get_from_*` family of static factory methods
            tlst7 = TreeList.get_from_stream(open('treefile.tre', 'rU'), "newick")
            tlst8 = TreeList.get_from_path('sometrees.nexus', "nexus")
            tlst9 = TreeList.get_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")

            # can also call `read()` on a TreeList object; each read adds the
            # tree(s) found to the TreeList
            tlst10 = TreeList()
            tlst10.read(open('boot1.tre', 'rU'), "newick")
            tlst10.read_from_stream(open('boot2.tre', 'rU'), "newick") # same as above
            tlst10.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick")
            tlst10.read_from_path("boot3.tre", "newick")

        """
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        iosys.Readable.__init__(self)
        iosys.Writeable.__init__(self)
        list.__init__(self)
        if len(args) > 1:
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        elif len(args) == 1:
            if hasattr(args[0], "__iter__") and not isinstance(args[0], str):
                stream = kwargs.get("stream")
                schema = kwargs.get("schema")
                if (stream is not None) or (schema is not None):
                    raise error.MultipleInitializationSourceError(class_name=self.__class__.__name__, arg=args[0])
                if isinstance(args[0], TreeList):
                    for t in args[0]:
                        if not isinstance(t, Tree):
                            raise ValueError("TreeList() requires TreeList or list of Tree objects as initialization argument in this context")
                        self.append(Tree(t))
                else:
                    for t in args[0]:
                        if not isinstance(t, Tree):
                            raise ValueError("TreeList() requires TreeList or list of Tree objects as initialization argument in this context (iterable including elements of type '%s' found instead)" % type(t))
                        self.append(t)
            else:
                raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        else:
            self.process_source_kwargs(**kwargs)

        if "oid" in kwargs:
            self.oid = kwargs["oid"]
        if "label" in kwargs:
            self.label = kwargs["label"]

    def __deepcopy__(self, memo):
        # we treat the taxa as immutable and copy the reference even in a deepcopy
        o = TaxonSetLinked.__deepcopy__(self, memo)
        for k, v in self.__dict__.iteritems():
            if k not in ['taxon_set', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

    def read(self, stream, schema, **kwargs):
        """
        Populates the `TreeList` from a `schema`-formatted file-like source
        `stream`. `schema` must be a recognized and tree file schema, such as
        `nexus`, `newick`, etc, for which a reader is available. If this is
        not implemented for the schema specified, then a
        `UnsupportedSchemaError` is raised. If the source defines multiple
        tree collections (e.g. multiple NEXUS "Trees" blocks), then the
        keyword argument ``collection_offset`` can be used to specify the
        0-based index of the tree collection. If not specified, or if value <
        0, then all collections in the source are merged and read.

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
        of the reader specialized to handle `schema` formats.
        """
        from dendropy.dataobject.dataset import DataSet
        collection_offset = kwargs.get("collection_offset", -1)
        tree_offset = kwargs.get("tree_offset", 0)
        if "taxon_set" in kwargs:
            if kwargs["taxon_set"] is not self.taxon_set and len(self) > 0:
                raise Exception("Cannot specify a different TaxonSet when reading into a populated TreeList.")
            else:
                self.taxon_set = kwargs["taxon_set"]
        else:
            kwargs["taxon_set"] = self.taxon_set
        kwargs["exclude_chars"] = True
        kwargs["exclude_trees"] = False
        d = DataSet(stream=stream, schema=schema, **kwargs)
        if len(d.tree_lists) == 0:
            raise ValueError("No trees in data source")
        if collection_offset >= len(d.tree_lists):
            raise IndexError("Tree collection offset %d specified, but data source only has %d tree collections defined" \
                % (collection_offset, len(d.tree_lists)))
        if collection_offset < 0:
            i = 0
            if self.label is None and len(self) == 0 and len(d.tree_lists) == 1:
                self.label = d.tree_lists[0].label
            for tlist in d.tree_lists:
                for t in tlist:
                    if i >= tree_offset:
                        self.append(t, reindex_taxa=False)
                    i += 1
            if i < tree_offset:
                raise IndexError("Tree offset %d specified, but data source only has %d trees defined" \
                        % (tree_offset, len(tlist)))
        else:
            tlist = d.tree_lists[collection_offset]
            if self.label is None and len(self) == 0 and tlist.label is not None:
                self.label = tlist.label
            if tree_offset < len(tlist):
                for t in tlist[tree_offset:]:
                    self.append(t, reindex_taxa=False)
            else:
                raise IndexError("Tree offset %d specified, but tree collection only has %d trees defined" \
                        % (tree_offset, len(tlist)))

    def write(self, stream, schema, **kwargs):
        """
        Writes out `TreeList` in `schema` to a destination described by
        given by `stream`.

        Additionally, for some formats, the following keywords are recognized:

            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
        """
        from dendropy.dataobject.dataset import DataSet
        d = DataSet()
        d.add(self)
        d.write(stream=stream, schema=schema, **kwargs)

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

    def consensus(self, min_freq=0.5, trees_splits_encoded=False, **kwargs):
        """
        Returns a consensus tree of all trees in self, with minumum frequency
        of split to be added to the consensus tree given by `min_freq`.
        """
        from dendropy import treesum
        self.split_distribution = treesplit.SplitDistribution(taxon_set=self.taxon_set)
        tsum = treesum.TreeSummarizer(**kwargs)
        tsum.count_splits_on_trees(self,
                split_distribution=self.split_distribution,
                trees_splits_encoded=trees_splits_encoded)
        tree = tsum.tree_from_splits(self.split_distribution, min_freq=min_freq)
        return tree

    def frequency_of_split(self, **kwargs):
        """
        Given a split or bipartition specified as:

            - a split bitmask given the keyword 'split_bitmask'
            - a list of `Taxon` objects given with the keyword `taxa`
            - a list of taxon labels given with the keyword `labels`
            - a list of oids given with the keyword `oids`

        this function returns the proportion of trees in self
        in which the split is found.
        """
        if "split_bitmask" in kwargs:
            split = kwargs["split_bitmask"]
        else:
            split = self.taxon_set.get_taxa_bitmask(**kwargs)
            k = kwargs.values()[0]
            if treesplit.count_bits(split) != len(k):
                raise IndexError('Not all taxa could be mapped to split (%s): %s' \
                    % (self.taxon_set.split_bitmask_string(split), k))
        found = 0
        total = 0
        for tree in self:
            if not hasattr(tree, "split_edges"):
                treesplit.encode_splits(tree)
            total += 1
            if split in tree.split_edges:
                found += 1
        return float(found)/total

    def __str__(self):
        return " ".join([ (str(tree) + ";") for tree in self ])

    def __repr__(self):
        return "<TreeList object at %s>" % (hex(id(self)))

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sTreeList object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            output_strio.write(':  %d Trees' % len(self))
            if depth >= 2:
                if self.taxon_set is not None:
                    tlead = "\n%s[Taxon Set]\n" % (" " * (indent+4))
                    output_strio.write(tlead)
                    self.taxon_set.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                tlead = "\n%s[Trees]\n" % (" " * (indent+4))
                output_strio.write(tlead)
                for i, t in enumerate(self):
                    t.description(depth=depth-1, indent=indent+8, itemize="[%d] " % (i), output=output_strio)
                    output_strio.write('\n')
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def as_python_source(self, tree_list_name=None, tree_list_args=None, oids=False):
        """
        Returns string that will rebuild this tree list in Python.
        """
        p = []

        if tree_list_name is None:
            tree_list_name = "tree_list_%s" % id(self)


        if self.label is not None:
            label = "'" + self.label + "'"
        else:
            label = "None"
        if oids:
            oid_str = ', oid="%s"' % self.oid
        else:
            oid_str = ""
        if tree_list_args is None:
            tree_list_args = ""
        else:
            tree_list_args = ", " + tree_list_args
        p.append("%s = dendropy.TreeList(label=%s%s%s)" \
            % (tree_list_name,
               label,
               oid_str,
               tree_list_args))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        taxon_map = {}
        for taxon in self.taxon_set:
            tobj_name = taxon_obj_namer(taxon)
            if taxon.label is not None:
                label = "'" + taxon.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % taxon.oid
            else:
                oid_str = ""
            p.append("%s = %s.taxon_set.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_list_name,
                   label,
                   oid_str))
            taxon_map[taxon] = tobj_name

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for tree in self:
            tree_obj_name = "tree_%s" % id(tree)
            if tree.label is not None:
                label = "'" + tree.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % tree.oid
            else:
                oid_str = ""
            p.append("%s = dendropy.Tree(label=%s, taxon_set=%s.taxon_set%s)" \
                % (tree_obj_name,
                   label,
                   tree_list_name,
                   oid_str))
            p.append("%s.append(%s, reindex_taxa=False)" % (tree_list_name, tree_obj_name))
            if oids:
                p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, tree.seed_node.oid))
            for node in tree.preorder_node_iter():
                for child in node.child_nodes():
                    if node is tree.seed_node:
                        nn = "%s.seed_node" % tree_obj_name
                    else:
                        nn = node_obj_namer(node)
                    if child.label is not None:
                        label = "'" + child.label + "'"
                    else:
                        label = "None"
                    if child.taxon is not None:
                        ct = taxon_obj_namer(child.taxon)
                    else:
                        ct = "None"
                    if oids:
                        oid_str = ', oid="%s"' % child.oid
                    else:
                        oid_str = ""
                    p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                            (node_obj_namer(child),
                             nn,
                             label,
                             ct,
                             child.edge.length,
                             oid_str))
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
    ## Special/Lifecycle methods

    def __init__(self, *args, **kwargs):
        """
        __init__ creates a new Tree object, optionally constructing it by cloning
        another Tree object if this is passed as the first argument, or
        out of a data source if `stream` and `schema` are keyword arguments are
        passed with a file-like object and a schema-specification string object
        values respectively.

        If `stream` and `schema` keyword arguments are given, will
        construct this `Tree` object from `schema`-formatted source
        given by file-like object `stream`. `schema` must be a
        recognized and tree file schema, such as `nexus`, `newick`, etc,
        for which a specialized tree list writer is available. If this
        is not implemented for the schema specified, then a
        `UnsupportedSchemaError` is raised. Other keywords will be
        passed to the underlying tree parser.

        Tree objects can thus be instantiated in the following ways::

            # /usr/bin/env python

            from cStringIO import StringIO
            from dendropy import Tree, TaxonSet

            # empty tree
            t1 = Tree()

            # tree from data source
            t2 = dendropy.Tree(stream=StringIO("((A,B),(C,D));"), schema="newick")

            # passing keywords to underlying tree parser
            t3 = dendropy.Tree(stream=StringIO("((A,B),(C,D));"),
                      schema="newick",
                      taxon_set=t3.taxon_set,
                      encode_splits=True)

            # tree structure deep-copied from another tree
            t4 = dendropy.Tree(t3)
            assert t4 is not t3                             # Trees are distinct
            assert t4.symmetric_difference(t3) == 0         # and structure is identical
            assert t4.taxon_set is t3.taxon_set             # BUT taxa are not cloned.
            nds3 = [nd for nd in t3.postorder_node_iter()]  # Nodes in the two trees
            nds4 = [nd for nd in t4.postorder_node_iter()]  # are distinct objects,
            for i, n in enumerate(nds3):                    # and can be manipulated
                assert nds3[i] is not nds4[i]               # independentally.
            egs3 = [eg for eg in t3.postorder_edge_iter()]  # Edges in the two trees
            egs4 = [eg for eg in t4.postorder_edge_iter()]  # are also distinct objects,
            for i, e in enumerate(egs3):                    # and can also be manipulated
                assert egs3[i] is not egs4[i]               # independentally.
            lves3 = t3.leaf_nodes()                         # Leaf nodes in the two trees
            lves4 = t4.leaf_nodes()                         # are also distinct objects,
            for i, lf in enumerate(lves3):                  # but order is the same,
                assert lves3[i] is not lves4[i]             # and associated Taxon objects
                assert lves3[i].taxon is lves4[i].taxon     # are the same.

            # to create deep copy of a tree with a different taxon set
            taxa = TaxonSet()
            t5 = dendropy.Tree(t3, taxon_set=taxa)
            assert t5 is not t3                             # As above, the trees are distinct
            assert t5.symmetric_difference(t3) == 0         # and the structures are identical,
            assert t5.taxon_set is not t3.taxon_set         # but this time, the taxa *are* different
            assert t5.taxon_set is taxa                     # as the given TaxonSet is used instead.
            lves3 = t3.leaf_nodes()                         # Leaf nodes (and, for that matter other nodes
            lves5 = t5.leaf_nodes()                         # as well as edges) are also distinct objects
            for i, lf in enumerate(lves3):                  # and the order is the same, as above,
                assert lves3[i] is not lves5[i]             # but this time the associated Taxon
                assert lves3[i].taxon is not lves5[i].taxon # objects are distinct though the taxon
                assert lves3[i].taxon.label == lves5[i].taxon.label # labels are the same.

            # the canonical way to instantiate a Tree from a data source
            # is the use the `get_from_*` family of static factory methods
            t6 = Tree.get_from_stream(open('treefile.tre', 'rU'), "newick", tree_offset=0)
            t7 = Tree.get_from_path('sometrees.nexus',
                    "nexus",
                    collection_offset=2,
                    tree_offset=1)
            s = "((A,B),(C,D));((A,C),(B,D));"
            t8 = Tree.get_from_string(s, "newick") # tree will be '((A,B),(C,D))'
            t9 = Tree.get_from_string(s, "newick", tree_offset=1) # tree will be '((A,C),(B,D))'

            # can also call `read()` on a Tree object; each read adds the
            # *replaces* the current tree with the definition specified in the
            # data source
            t10 = Tree()
            t10.read(open('boot1.tre', 'rU'), "newick", tree_offset=0)
            t10.read_from_stream(open('boot2.tre', 'rU'), "newick") # same as above
            t10.read_from_string("((A,B),(C,D));((A,C),(B,D));", "newick", tree_offset=0)
            t10.read_from_path("mle.tre", "newick")

            # to 'switch out' the TaxonSet of a tree, replace the reference and
            # reindex the taxa:
            t11 = Tree.get_from_string('((A,B),(C,D));', 'newick')
            taxa = TaxonSet()
            t11.taxon_set = taxa
            t11.reindex_subcomponent_taxa()

        """
        TaxonSetLinked.__init__(self,
                                taxon_set=kwargs.get("taxon_set", None),
                                label=kwargs.get("label", None),
                                oid=kwargs.get("oid", None))
        iosys.Writeable.__init__(self)
        iosys.Readable.__init__(self)
        self.seed_node = Node(edge=Edge())
        self.length_type = None
        self.comments = None
        self._is_rooted = None
        self.weight = None

        if len(args) > 1:
            raise error.TooManyArgumentsError(func_name=self.__class__.__name__, max_args=1, args=args)
        if len(args) == 1:
            if ("stream" in kwargs and kwargs["stream"] is not None) \
                    or ("schema" in kwargs and kwargs["schema"] is not None):
                raise error.MultipleInitializationSourceError(class_name=self.__class__.__name__, arg=args[0])
            if isinstance(args[0], Node):
                self.seed_node = args[0]
            elif isinstance(args[0], Tree):
                self.clone_from(args[0])
                if "taxon_set" in kwargs:
                    self.taxon_set = kwargs["taxon_set"]
                    self.reindex_subcomponent_taxa()
            else:
                raise error.InvalidArgumentValueError(func_name=self.__class__.__name__, arg=args[0])
        else:
            if "seed_node" in kwargs:
                self.seed_node = kwargs['seed_node']
            self.process_source_kwargs(**kwargs)
        if "oid" in kwargs:
            self.oid = kwargs["oid"]
        if "label" in kwargs:
            self.label = kwargs["label"]

    ###########################################################################
    ## I/O

    def clone_from(self, other):
        """
        Clones the structure and properties of `Tree` object `other`.
        """
        t = copy.deepcopy(other)
        self.__dict__ = t.__dict__
        return self

    def __deepcopy__(self, memo):
        # we treat the taxa as immutable and copy the reference even in a deepcopy
        o = TaxonSetLinked.__deepcopy__(self, memo)
        for k, v in self.__dict__.iteritems():
            if k not in ['taxon_set', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        return o

    def read(self, stream, schema, **kwargs):
        """
        Populates/constructs objects of this type from `schema`-formatted
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

        If the source defines multiple tree collections (e.g. multiple NEXUS "Trees"
        blocks), then the keyword argument ``collection_offset`` can be used to specify
        the 0-based index of the tree collection, and the keyword argument
        ``tree_offset`` can be used to specify the 0-based index of the tree
        within the collection, as the source.
        If ``collection_offset`` is not specified and < 0, then all collections
        in the source are merged before considering ``tree_offset``.
        If ``tree_offset`` is not specified, then the first tree (offset=0) is
        returned.
        """
        from dendropy.dataobject.dataset import DataSet
        collection_offset = kwargs.get("collection_offset", -1)
        tree_offset = kwargs.get("tree_offset", 0)
        if "taxon_set" in kwargs:
            self.taxon_set = kwargs["taxon_set"]
        kwargs["exclude_chars"] = True
        kwargs["exclude_trees"] = False
        d = DataSet(stream=stream, schema=schema, **kwargs)
        if len(d.tree_lists) == 0:
            raise ValueError("No trees in data source")
        if collection_offset >= len(d.tree_lists):
            raise IndexError("Tree collection offset %d specified, but data source only has %d tree collections defined" \
                % (collection_offset, len(d.tree_lists)))
        if collection_offset < 0:
            i = 0
            for tlist in d.tree_lists:
                for t in tlist:
                    if i == tree_offset:
                        self.__dict__ = t.__dict__
                    i += 1
            if i <= tree_offset:
                raise IndexError("Tree offset %d specified, but data source only has %d trees defined" \
                            % (tree_offset, i+1))
        else:
            tlist = d.tree_lists[collection_offset]
            if tree_offset < len(tlist):
                t = tlist[tree_offset]
                self.__dict__ = t.__dict__
            else:
                raise IndexError("Tree offset %d specified, but tree collection only has %d trees defined" \
                        % (tree_offset, len(tlist)))

    def write(self, stream, schema, **kwargs):
        """
        Writes out `Tree` in `schema` to a destination given by file-like object
        `stream`.

        `schema` must be a recognized and tree file schema, such as `nexus`,
        `newick`, etc, for which a specialized tree list writer is
        available. If this is not implemented for the schema specified, then
        a `UnsupportedSchemaError` is raised.

        Additionally, for some formats, the following keywords are recognized:

            - `edge_lengths` : if False, edges will not write edge lengths [True]
            - `internal_labels` : if False, internal labels will not be written [True]
        """
        from dendropy.dataobject.dataset import DataSet
        d = DataSet()
        d.add(TreeList([self], taxon_set=self.taxon_set))
        d.write(stream=stream, schema=schema, **kwargs)

    ###########################################################################
    ## Accessors

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

    def find_node_with_label(self, label):
        "Finds the first node with matching label."
        for node in self.preorder_node_iter():
            if node.label == label:
                return node
        return None

    def find_node_with_taxon(self, taxon_filter_fn=None):
        "Finds the first node for which taxon_filter_fn(node.taxon) == True."
        for node in self.preorder_node_iter():
            if hasattr(node, "taxon") and node.taxon is not None:
                if taxon_filter_fn(node.taxon):
                    return node
        return None

    def find_node_with_taxon_label(self, label):
        "Returns node with taxon with given label."
        taxon = self.taxon_set.get_taxon(label=label)
        if taxon is None:
            return None
        return self.find_node_with_taxon(lambda x: x is taxon)

    def find_edge(self, oid):
        "Finds the first edge with matching id."
        for edge in self.preorder_edge_iter():
            if edge.oid == oid:
                return edge
        return None

    def get_edge_set(self, filter_fn=None):
        """
        Returns the set of edges that are currently in the tree.
        Note: the returned set acts like a shallow copy of the edge set (adding
        or deleting elements from the set does not change the tree, but
        modifying the elements does).
        """
        return set([i for i in self.preorder_edge_iter(filter_fn=filter_fn)])

    def get_node_set(self, filter_fn=None):
        """Returns the set of nodes that are currently in the tree

        Note: the returned set acts like a shallow copy of the edge set (adding
        or deleting elements from the set does not change the tree, but
        modifying the elements does).
        """
        return set([i for i in self.preorder_node_iter(filter_fn=filter_fn)])

    def mrca(self, **kwargs):
        """
        Returns the shallowest node in the tree (the node furthest from
        the root, or `start_node`, in the direction toward the tips of
        the tree) that has all of the taxa that:

            - are specified by the split bitmask given by the keyword argument `split_bitmask`
            - are in the list of Taxon objects given by the keyword argument 'taxa'
            - have the labels specified by the list of strings given by the keyword argument 'taxon_labels'

        Returns None if no appropriate node is found.
        Assumes that edges on tree have been decorated with treesplit.
        It is possible that split is not compatible with the subtree that is
        returned! (compatibility tests are not fully performed).
        This function is used to find the "insertion point" for a new split via a
        root to tip search.
        """
        start_node = kwargs.get("start_node", self.seed_node)
        split_bitmask = None
        if "split_bitmask" in kwargs:
            split_bitmask = kwargs["split_bitmask"]
        else:
            taxa = kwargs.get("taxa", None)
            if taxa is None:
                if "taxon_labels" in kwargs:
                    taxa = self.taxon_set.get_taxa(labels=kwargs["taxon_labels"])
                    if len(taxa) != len(kwargs["taxon_labels"]):
                        raise KeyError("Not all labels matched to taxa")
                else:
                    raise TypeError("Must specify one of: 'split_bitmask', 'taxa' or 'taxon_labels'")
            if taxa is None:
                raise ValueError("No taxa matching criteria found")
            split_bitmask = self.taxon_set.get_taxa_bitmask(taxa=taxa)

        if split_bitmask is None or split_bitmask == 0:
            raise ValueError("Null split bitmask (0)")

        if not hasattr(start_node.edge, "split_bitmask"):
            treesplit.encode_splits(self, delete_outdegree_one=False)

        if (start_node.edge.split_bitmask & split_bitmask) != split_bitmask:
            return None

        curr_node = start_node
        last_match = start_node
        nd_source = iter(start_node.child_nodes())
        try:
            while True:
                cm = curr_node.edge.split_bitmask
                cms = (cm & split_bitmask)
                if cms:
                    # for at least one taxon cm has 1 and split has 1
                    if cms == split_bitmask:
                        # curr_node has all of the 1's that split has
                        if cm == split_bitmask:
                            return curr_node
                        last_match = curr_node
                        nd_source = iter(curr_node.child_nodes())
                    else:
                        # we have reached a child that has some, but not all of the
                        #   required taxa as descendants, so we return the last_match
                        return last_match
                curr_node = nd_source.next()
        except StopIteration:
            # we shouldn't reach this if all of the descendants are properly
            #   decorated with split_bitmask attributes, but there may be some hacky
            #   context in which we want to allow the function to be called with
            #   leaves that have not been encoded with split_bitmasks.
            return last_match

    ###########################################################################
    ## Node iterators

    def __iter__(self):
        return self.preorder_node_iter()

    def preorder_node_iter(self, filter_fn=None):
        "Returns preorder iterator over tree nodes."
        for node in self.seed_node.preorder_iter(filter_fn=filter_fn):
            yield node

    def postorder_node_iter(self, filter_fn=None):
        "Returns postorder iterator over tree nodes."
        for node in self.seed_node.postorder_iter(filter_fn=filter_fn):
            yield node

    def level_order_node_iter(self, filter_fn=None):
        "Returns level-order iterator over tree nodes."
        for node in self.seed_node.level_order_iter(filter_fn=filter_fn):
            yield node

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over tree leaf_nodes (order determined by
        postorder tree-traversal).
        """
        for node in self.seed_node.leaf_iter(filter_fn=filter_fn):
            yield node

    def age_order_node_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Iterates over nodes in order of age. If `include_leaves` is False, will
        skip leaves (default is not to skip leaves). If `descending` is True,
        will go from oldest nodes to youngest (default is asecending: youngest
        nodes to oldest).
        """
        if self.seed_node.age is None:
            self.calc_node_ages()
        for node in self.seed_node.age_order_iter(include_leaves=include_leaves, filter_fn=filter_fn, descending=descending):
            yield node

    def postorder_internal_node_iter(self, filter_fn=None):
        """
        Iterates over all internal nodes in post-order.
        """
        for node in self.seed_node.postorder_internal_node_iter(filter_fn=filter_fn):
            yield node

    def preorder_internal_node_iter(self, filter_fn=None):
        """
        Iterates over all internal nodes in pre-order.
        """
        for node in self.seed_node.preorder_internal_node_iter(filter_fn=filter_fn):
            yield node

    ###########################################################################
    ## Edge iterators

    def preorder_edge_iter(self, filter_fn=None):
        "Returns preorder iterator over tree edges."
        for node in self.seed_node.preorder_iter():
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def postorder_edge_iter(self, filter_fn=None):
        "Returns postorder iterator over tree edges."
        for node in self.seed_node.postorder_iter():
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def level_order_edge_iter(self, filter_fn=None):
        "Returns level-order iterator over tree edges."
        for node in self.seed_node.level_order_iter():
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    def leaf_edge_iter(self, filter_fn=None):
        "Returns iterator over tree leaf edges."
        for node in self.seed_node.leaf_iter():
            if node.edge and (filter_fn is None or filter_fn(node.edge)):
                yield node.edge

    ###########################################################################
    ## Taxa Management

    def infer_taxa(self):
        """
        Returns a new TaxonSet object populated with taxa from this
        tree.
        """
        taxon_set = TaxonSet()
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

    def unassign_taxa(self, exclude_leaves=False, exclude_internal=False):
        """
        Strips taxon assignments from tree. If `exclude_leaves` is True,
        then taxa on leaves will be retained. If `exclude_internal` is True,
        then taxa on internal nodes will be retained. The `taxon_set` is not
        affected by this operation.
        """
        for nd in self.postorder_node_iter():
            if (len(nd._child_nodes) == 0) and not exclude_leaves:
                nd.taxon = None
            elif (len(nd._child_nodes) > 0) and not exclude_internal:
                nd.taxon = None

    def randomly_assign_taxa(self, create_required_taxa=True, rng=None):
        """
        Randomly assigns taxa to leaf nodes. If the number of taxa defined in
        the taxon set of the tree is more than the number of tips, then a random
        subset of taxa in `taxon_set` will be assigned to the tips of tree.
        If the number of tips is more than the number of taxa in the `taxon_set`,
        and `add_extra_taxa` is not True [default], then new Taxon
        objects will be created and added to the `taxon_set`; if `create_required_taxa`
        is False, then an exception is raised.

        In addition, a Random() object or equivalent can be passed using `rng`;
        otherwise GLOBAL_RNG is used.
        """
        if rng is None:
            rng = GLOBAL_RNG
        if len(self.taxon_set) == 0:
            for i, nd in enumerate(self.leaf_nodes()):
                nd.taxon = self.taxon_set.require_taxon(label=("T%d" % (i+1)))
        else:
            taxa = [t for t in self.taxon_set]
            for i, nd in enumerate(self.leaf_nodes()):
                if len(taxa) > 0:
                    nd.taxon = taxa.pop(rng.randint(0, len(taxa)-1))
                else:
                    if not create_required_taxa:
                        raise ValueError("TaxonSet has %d taxa, but tree has %d tips" % (len(self.taxon_set), len(self.leaf_nodes())))
                    label = "T%d" % (i+1)
                    k = 0
                    while self.taxon_set.has_taxon(label=label):
                        label = "T%d" % (i+1+k)
                        k += 1
                    nd.taxon = self.taxon_set.require_taxon(label=label)

    ###########################################################################
    ## Structure

    def _get_rooting_state_is_undefined(self):
        return self._is_rooted is None

    rooting_state_is_undefined = property(_get_rooting_state_is_undefined)

    def _get_is_rooted(self):
        return self._is_rooted

    def _set_is_rooted(self, val):
        self._is_rooted = val

    is_rooted = property(_get_is_rooted, _set_is_rooted)

    def _get_is_unrooted(self):
        return not self._is_rooted

    def _set_is_unrooted(self, val):
        self._is_rooted = not val

    is_unrooted = property(_get_is_unrooted, _set_is_unrooted)

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
        to_del_edge.collapse()
        self.is_rooted = False
        return self.seed_node

    def reseed_at(self, new_seed_node, update_splits=False, delete_outdegree_one=True):
        """
        Takes an internal node, `new_seed_node` that must already be in the tree and
        rotates the tree such that `new_seed_node` is the `seed_node` of the tree.
        This is a 'soft' rerooting -- i.e., changes the tree representation so
        tree traversal behaves as if the tree is rooted at 'new_seed_node', but
        it does not actually change the tree's rooting state.  If
        `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        if new_seed_node.is_leaf():
            raise ValueError('Rooting at a leaf is not supported')

        old_par = new_seed_node.parent_node
        if old_par is None:
            return
        full_encode = False
        if update_splits:
            try:
                taxa_mask = self.seed_node.edge.split_bitmask
            except:
                full_encode = True
                update_splits = False
        to_edge_dict = None
        if update_splits:
            to_edge_dict = getattr(self, "split_edges", None)

        if old_par is self.seed_node:
            root_children = old_par.child_nodes()
            if len(root_children) == 2 and delete_outdegree_one:
                # root (old_par) was of degree 2, thus we need to suppress the
                #   node
                fc = root_children[0]
                if fc is new_seed_node:
                    sister = root_children[1]
                else:
                    assert root_children[1] is new_seed_node
                    sister = fc
                if new_seed_node.edge.length:
                    sister.edge.length += new_seed_node.edge.length
                edge_to_del = new_seed_node.edge
                new_seed_node.edge = old_par.edge
                if update_splits:
                    assert new_seed_node.edge.split_bitmask == taxa_mask
                if to_edge_dict:
                    del to_edge_dict[edge_to_del.split_bitmask]
                new_seed_node.add_child(sister, edge_length=sister.edge.length)
                self.seed_node = new_seed_node
                return
        else:
            self.reseed_at(old_par,
                    update_splits=update_splits,
                    delete_outdegree_one=delete_outdegree_one)
        old_par.edge, new_seed_node.edge = new_seed_node.edge, old_par.edge
        e = old_par.edge
        if update_splits:
            if to_edge_dict:
                del to_edge_dict[e.split_bitmask]
            e.split_bitmask = (~(e.split_bitmask)) & taxa_mask
            if to_edge_dict:
                to_edge_dict[e.split_bitmask] = e
            assert new_seed_node.edge.split_bitmask == taxa_mask
        old_par.remove_child(new_seed_node)
        new_seed_node.add_child(old_par, edge_length=e.length)
        self.seed_node = new_seed_node
        if full_encode:
            treesplit.encode_splits(self, delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def to_outgroup_position(self, outgroup_node, update_splits=False, delete_outdegree_one=True):
        """Reroots the tree at the parent of `outgroup_node` and makes `outgroup_node` the first child
        of the new root.  This is just a convenience function to make it easy
        to place a clade as the first child under the root.
        Assumes that `outgroup_node` and `outgroup_node.parent_node` and are in the tree/
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        p = outgroup_node.parent_node
        assert p is not None
        self.reseed_at(p, update_splits=update_splits, delete_outdegree_one=delete_outdegree_one)
        p.remove_child(outgroup_node)
        p.add_child(outgroup_node, edge_length=outgroup_node.edge.length, pos=0)
        return self.seed_node

    def reroot_at_node(self, new_root_node, update_splits=False, delete_outdegree_one=True):
        """
        Takes an internal node, `new_seed_node` that must already be in the tree and
        roots the tree at that node.
        This is a 'hard' rerooting -- i.e., changes the tree
        representation so tree traversal behaves as if the tree is rooted at
        'new_seed_node', *and* changes the tree's rooting state.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        self.reseed_at(new_seed_node=new_root_node,
                update_splits=False,
                delete_outdegree_one=delete_outdegree_one)
        self.is_rooted = True
        if update_splits:
            self.update_splits(delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def reroot_at_edge(self,
            edge,
            length1=None,
            length2=None,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Takes an internal edge, `edge`, adds a new node to it, and then roots
        the tree on the new node.
        `length1` and `length2` will be assigned to the new (sub-)edge leading
        to the old parent of the original edge, while `length2` will be
        assigned to the old child of the original edge.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        old_tail = edge.tail_node
        old_head = edge.head_node
        new_seed_node = old_tail.new_child(edge_length=length1)
        old_tail.remove_child(old_head)
        new_seed_node.add_child(old_head, edge_length=length2)
        self.reroot_at_node(new_seed_node,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)
        return self.seed_node

    def reroot_at_midpoint(self, update_splits=False, delete_outdegree_one=True):
        """
        Reroots the tree at the the mid-point of the longest distance between
        two taxa in a tree.
        Sets the rooted flag on the tree to True.
        If `update_splits` is True, then the edges' `split_bitmask` and the tree's
        `split_edges` attributes will be updated.
        If the *old* root of the tree had an outdegree of 2, then after this
        operation, it will have an outdegree of one. In this case, unless
        `delete_outdegree_one` is False, then it will be
        removed from the tree.
        """
        from dendropy import treecalc
        pdm = treecalc.PatristicDistanceMatrix(self)
        n1, n2 = pdm.max_dist_nodes
        plen = float(pdm.max_dist) / 2
        mrca_node = pdm.mrca(n1.taxon, n2.taxon)
        #assert mrca_node is self.mrca(taxa=[n1.taxon, n2.taxon])
        #mrca_node = self.mrca(taxa=[n1.taxon, n2.taxon])
        cur_node = n1

        break_on_node = None # populated *iff* midpoint is exactly at an existing node
        target_edge = None
        head_node_edge_len = None

        # going up ...
        while cur_node is not mrca_node:
            if cur_node.edge.length > plen:
                target_edge = cur_node.edge
                head_node_edge_len = plen #cur_node.edge.length - plen
                plen = 0
                break
            elif cur_node.edge.length < plen:
                plen -= cur_node.edge.length
                cur_node = cur_node.parent_node
            else:
                break_on_node = cur_node

        assert break_on_node is not None or target_edge is not None

        if break_on_node:
            self.reseed_at(break_on_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
            new_seed_node = break_on_node
        else:
            tail_node_edge_len = target_edge.length - head_node_edge_len
            old_head_node = target_edge.head_node
            old_tail_node = target_edge.tail_node
            old_tail_node.remove_child(old_head_node)
            new_seed_node = Node()
            new_seed_node.add_child(old_head_node, edge_length=head_node_edge_len)
            old_tail_node.add_child(new_seed_node, edge_length=tail_node_edge_len)
            self.reseed_at(new_seed_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
        self.is_rooted = True
        if update_splits:
            self.update_splits(delete_outdegree_one=False)
        return self.seed_node

    def delete_outdegree_one_nodes(self):
        for nd in self.postorder_node_iter():
            children = nd.child_nodes()
            if len(children) == 1:
                if nd.edge.length is not None:
                    if children[0].edge.length is None:
                        children[0].edge.length = nd.edge.length
                    else:
                        children[0].edge.length += nd.edge.length
                if nd.parent_node is not None:
                    pos = nd.parent_node.child_nodes().index(nd)
                    nd.parent_node.add_child(children[0], pos=pos)
                    nd.parent_node.remove_child(nd)
                else:
                    assert nd is self.seed_node
                    self.seed_node = children[0]
                    self.seed_node.parent_node = None

    def collapse_unweighted_edges(self,
            threshold=0.0000001,
            update_splits=False):
        """
        Collapse all edges with edge lengths less than or equal to
        ``threshold``.
        """
        for e in self.postorder_edge_iter():
            if e.length <= threshold:
               e.collapse()
        if update_splits:
            self.update_splits()

    def resolve_polytomies(self, update_splits=False, rng=None):
        """
        Arbitrarily resolve polytomies using 0-length splits.

        If `rng` is an object with a sample() method then the polytomy will be
            resolved by sequentially adding (generating all tree topologies
            equiprobably
            rng.sample() should behave like random.sample()
        If `rng` is not passed in, then polytomy is broken deterministically by
            repeatedly joining pairs of children.
        """
        polytomies = []
        for node in self.postorder_node_iter():
            if len(node.child_nodes()) > 2:
                polytomies.append(node)
        for node in polytomies:
            children = node.child_nodes()
            nc = len(children)
            if nc > 2:
                if rng:
                    to_attach = children[2:]
                    for child in to_attach:
                        node.remove_child(child)
                    attachment_points = children[:2] + [node]
                    while len(to_attach) > 0:
                        next_child = to_attach.pop()
                        next_sib = rng.sample(attachment_points, 1)[0]
                        next_attachment = Node()
                        p = next_sib.parent_node
                        p.add_child(next_attachment)
                        p.remove_child(next_sib)
                        next_attachment.add_child(next_sib)
                        next_attachment.add_child(next_child)
                        attachment_points.append(next_attachment)
                else:
                    while len(children) > 2:
                        nn1 = Node()
                        nn1.edge.length = 0
                        c1 = children[0]
                        c2 = children[1]
                        node.remove_child(c1)
                        node.remove_child(c2)
                        nn1.add_child(c1)
                        nn1.add_child(c2)
                        node.add_child(nn1)
                        children = node.child_nodes()
        if update_splits:
            self.update_splits()

    def prune_subtree(self,
            node,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes subtree starting at `node` from tree.
        """
        if not node:
            raise ValueError("Tried to remove an non-existing or null node")
        if node.parent_node is None:
            raise TypeError('Node has no parent and is implicit root: cannot be pruned')
        node.parent_node.remove_child(node)
        if delete_outdegree_one:
            self.delete_outdegree_one_nodes()
        if update_splits:
            self.update_splits()

    def prune_leaves_without_taxa(self,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes all terminal nodes that have their ``taxon`` attribute set to
        ``None``.
        """
        for nd in self.leaf_iter():
            if nd.taxon is None:
                nd.edge.tail_node.remove_child(nd)
        if delete_outdegree_one:
            self.delete_outdegree_one_nodes()
        if update_splits:
            self.update_splits()

    def prune_taxa(self, taxa, update_splits=False, delete_outdegree_one=True):
        """
        Removes terminal nodes associated with Taxon objects given by the container
        `taxa` (which can be any iterable, including a TaxonSet object) from `self`.
        """
        nodes = []
        for taxon in taxa:
            nd = self.find_node(lambda x: x.taxon is taxon)
            if nd is not None:
                nd.edge.tail_node.remove_child(nd)
        self.prune_leaves_without_taxa(update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def prune_taxa_with_labels(self,
            labels,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are associated with Taxon objects with
        labels given by `labels`.
        """
        taxa = self.taxon_set.get_taxa(labels=labels)
        self.prune_taxa(taxa=taxa,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def retain_taxa(self,
            taxa,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are not associated with any
        of the Taxon objects given by ``taxa`` (which can be any iterable, including a
        TaxonSet object) from the ``self``.
        """
        to_prune = [t for t in self.taxon_set if t not in taxa]
        self.prune_taxa(to_prune,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def retain_taxa_with_labels(self,
            labels,
            update_splits=False,
            delete_outdegree_one=True):
        """
        Removes terminal nodes that are not associated with Taxon objects with
        labels given by `labels`.
        """
        taxa = self.taxon_set.get_taxa(labels=labels)
        self.retain_taxa(taxa=taxa,
                update_splits=update_splits,
                delete_outdegree_one=delete_outdegree_one)

    def randomly_reorient_tree(self, rng=None, update_splits=False):
        """
        Randomly picks a new rooting position and rotates the branches around all
        internal nodes in the `self`. If `update_splits` is True, the the `split_bitmask`
        and `split_edges` attributes kept valid.
        """
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        nd = rng.sample(self.nodes(), 1)[0]
        if nd.is_leaf():
            self.to_outgroup_position(nd, update_splits=update_splits)
        else:
            self.reseed_at(nd, update_splits=update_splits)
        self.randomly_rotate(rng=rng)

    def randomly_rotate(self, rng=None):
        "Randomly rotates the branches around all internal nodes in `self`"
        if rng is None:
            rng = GLOBAL_RNG # use the global rng by default
        internal_nodes = self.internal_nodes()
        for nd in internal_nodes:
            c = nd.child_nodes()
            rng.shuffle(c)
            nd.set_child_nodes(c)

    def ladderize(self, ascending=True):
        """
        Sorts child nodes in ascending (if ``ascending`` is ``False``) or
        descending (if ``ascending`` is ``False``) order in terms of the number of
        children each child node has.
        """
        node_desc_counts = {}
        for nd in self.postorder_node_iter():
            if len(nd._child_nodes) == 0:
                node_desc_counts[nd] = 0
            else:
                total = 0
                for child in nd._child_nodes:
                    total += node_desc_counts[child]
                total += len(nd._child_nodes)
                node_desc_counts[nd] = total
                nd._child_nodes.sort(key=lambda n: node_desc_counts[n], reverse=not ascending)

    def update_splits(self, **kwargs):
        """
        Recalculates split hashes for tree.
        """
        treesplit.encode_splits(self, **kwargs)

    ###########################################################################
    ## Ages, depths, branch lengths etc. (mutation)

    def scale_edges(self, edge_len_multiplier):
        """Multiplies every edge length in `self` by `edge_len_multiplier`"""
        for e in self.postorder_edge_iter():
            if e.length is not None:
                e.length *= edge_len_multiplier

    def set_edge_lengths_from_node_ages(self, allow_negative_edges=False):
        """
        Sets the edge lengths of the tree so that the path lengths from the
        tips equal the value of the `age` attribute of the nodes.
        """
        for nd in self.preorder_node_iter():
            if nd.parent_node is not None:
                #if nd.parent_node.age < nd.age:
                #    nd.edge.length = 0.0
                #else:
                #    nd.edge.length = nd.parent_node.age - nd.age
                if not allow_negative_edges and nd.parent_node.age < nd.age:
                    #if nd.parent_node is self.seed_node:
                    #    # special case seed node
                    #    nd.parent_node.age = nd.age + nd.edge_length
                    #else:
                    #    raise ValueError('Parent node age (%s: %s) is younger than descendent (%s: %s)'
                    #            % (nd.parent_node.oid, nd.parent_node.age, nd.oid, nd.age))
                    raise ValueError('Parent node age (%s: %s) is younger than descendent (%s: %s)'
                            % (nd.parent_node.oid, nd.parent_node.age, nd.oid, nd.age))
                nd.edge.length = nd.parent_node.age - nd.age

    ###########################################################################
    ## Ages, depths, branch lengths etc. (calculation)

    def calc_node_ages(self, check_prec=0.0000001):
        """
        Adds an attribute called "age" to  each node, with the value equal to
        the sum of edge lengths from the node to the tips. If the lengths of
        different paths to the node differ by more than `check_prec`, then a
        ValueError exception will be raised indicating deviation from
        ultrametricity. If `check_prec` is negative or False, then this check
        will be skipped.
        """
        for node in self.postorder_node_iter():
            ch = node.child_nodes()
            if len(ch) == 0:
               node.age = 0.0
            else:
                first_child = ch[0]
                node.age = first_child.age + first_child.edge.length
                if not (check_prec < 0 or check_prec == False):
                    for nnd in ch[1:]:
                        ocnd = nnd.age + nnd.edge.length
                        if abs(node.age - ocnd) > check_prec:
                            raise ValueError("Tree is not ultrametric")

    def node_ages(self, check_prec=0.0000001):
        """
        Returns list of ages of speciation events / coalescence times on tree.
        """
        try:
            ages = [n.age for n in self.internal_nodes()]
        except AttributeError:
            self.calc_node_ages(check_prec=check_prec)
            ages = [n.age for n in self.internal_nodes()]
        ages.sort()
        return ages

    def length(self):
        """
        Returns sum of edge lengths of self. Edges with no lengths defined
        (None) will be considered to have a length of 0.
        Note that we do not overrride `__len__` as this requires an integer
        return value.
        """
        total = 0
        for edge in self.postorder_edge_iter():
            if edge.length is not None:
                total += edge.length
        return total

    def coalescence_intervals(self):
        """
        Returns list of coalescence intervals of self., i.e., the waiting
        times between successive coalescence events.
        """
        ages = self.node_ages()
        intervals = []
        intervals.append(ages[0])
        for i, d in enumerate(ages[1:]):
            intervals.append(d - ages[i])
        return intervals

    ###########################################################################
    ## Metrics -- Internal

    def pybus_harvey_gamma(self, prec=0.00001):
        """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
        is used to test for constancy of birth and death rates over the course of
        a phylogeny.  Under the pure-birth process, the statistic should follow
        a standard Normal distibution: a Normal(mean=0, variance=1).

        If the lengths of different paths to the node differ by more than `prec`,
            then a ValueError exception will be raised indicating deviation from
            ultrametricty.
        Raises a Value Error if the tree is not ultrametric, is non-binary, or has
            only 2 leaves.

        As a side effect a `age` attribute is added to the nodes of the self.

        Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
        molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
        (267). 2267-2272
        """
        # the equation is given by:
        #   T = \sum_{j=2}^n (jg_j)
        #   C = T \sqrt{\frac{1}{12(n-2)}}
        #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
        # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
        #   times between consecutive (in time, not along a branch) speciation times.
        node = None
        speciation_ages = []
        n = 0
        if self.seed_node.age is None:
            self.calc_node_ages(check_prec=prec)
        for node in self.postorder_node_iter():
            if len(node.child_nodes()) == 2:
                speciation_ages.append(node.age)
            else:
                n += 1
        if node is None:
            raise ValueError("Empty tree encountered")
        speciation_ages.sort(reverse=True)
        g = []
        older = speciation_ages[0]
        for age in speciation_ages[1:]:
            g.append(older - age)
            older = age
        g.append(older)
        if not g:
            raise ValueError("No internal nodes found (other than the root)")
        assert(len(g) == (n - 1))
        T = 0.0
        accum = 0.0
        for i in xrange(2, n):
            list_index = i - 2
            T += i * float(g[list_index])
            accum += T
        list_index = n - 2
        T += (n) * g[list_index]
        nmt = n - 2.0
        numerator = accum/nmt - T/2.0
        C = T*pow(1/(12*nmt), 0.5)
        return numerator/C

    ###########################################################################
    ## Metrics -- Comparative

    def find_missing_splits(self, other_tree):
        """
        Returns a list of splits that are in self,  but
        not in `other_tree`.
        """
        missing = []
        if self.taxon_set is not other_tree.taxon_set:
            raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                    % (hex(id(self.taxon_set)), hex(id(other_tree.taxon_set))))
        if not hasattr(self, "split_edges"):
            self.encode_splits()
        if not hasattr(other_tree, "split_edges"):
            other_tree.encode_splits()
        for split in self.split_edges:
            if split in other_tree.split_edges:
                pass
            else:
                missing.append(split)
        return missing

    def symmetric_difference(self, other_tree):
        """
        Returns the symmetric_distance between this tree and the tree given by
        `other`, i.e. the sum of splits found in one but not in both trees.
        """
        t = self.false_positives_and_negatives(other_tree)
        return t[0] + t[1]

    def false_positives_and_negatives(self, other_tree):
        """
        Returns a tuple pair: all splits found in `other` but in self, and all
        splits in self not found in other.
        """
        from dendropy import treecalc
        if other_tree.taxon_set is not self.taxon_set:
            other_tree = Tree(other_tree, taxon_set=self.taxon_set)
        return treecalc.false_positives_and_negatives(self, other_tree)

    def robinson_foulds_distance(self, other_tree):
        """
        Returns Robinson-Foulds distance between this tree and `other_tree`.
        """
        from dendropy import treecalc
        if other_tree.taxon_set is not self.taxon_set:
            other_tree = Tree(other_tree, taxon_set=self.taxon_set)
        return treecalc.robinson_foulds_distance(self, other_tree)

    def euclidean_distance(self, other_tree):
        """
        Returns Euclidean_distance distance between this tree and `other_tree`.
        """
        from dendropy import treecalc
        if other_tree.taxon_set is not self.taxon_set:
            other_tree = Tree(other_tree, taxon_set=self.taxon_set)
        return treecalc.euclidean_distance(self, other_tree)

    ###########################################################################
    ## Metadata

    def strip_comments(self):
        """
        Remove comments from tree/nodes.
        """
        self.comments = []
        for nd in self.postorder_node_iter():
            nd.comments = []
            nd.edge.comments = []

    ###########################################################################
    ## Representation

    def __str__(self):
        "Dump Newick string."
        return "%s" % self.as_newick_string()

    def __repr__(self):
        return "<Tree object at %s>" % (hex(id(self)))

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sTree object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            newick_str = self.as_newick_string()
            if not newick_str:
                newick_str = "()"
            if depth == 1:
                output_strio.write(': %s' % newick_str)
            elif depth >= 2:
                num_nodes = len([nd for nd in self.preorder_node_iter()])
                num_edges = len([ed for ed in self.preorder_edge_iter()])
                output_strio.write(': %d Nodes, %d Edges' % (num_nodes, num_edges))
                if self.taxon_set is not None:
                    output_strio.write("\n%s[Taxon Set]\n" % (" " * (indent+4)))
                    self.taxon_set.description(depth=depth-1, indent=indent+8, itemize="", output=output_strio)
                output_strio.write('\n%s[Tree]' % (" " * (indent+4)))
                output_strio.write('\n%s%s' % (" " * (indent+8), newick_str))
                if depth >= 3:
                    output_strio.write("\n%s[Nodes]" % (" " * (indent+4)))
                    for i, nd in enumerate(self.preorder_node_iter()):
                        output_strio.write('\n')
                        nd.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_set=self.taxon_set)
                    output_strio.write("\n%s[Edges]" % (" " * (indent+4)))
                    for i, ed in enumerate(self.preorder_edge_iter()):
                        output_strio.write('\n')
                        ed.description(depth=depth-3, indent=indent+8, itemize="[%d] " % i, output=output_strio, taxon_set=self.taxon_set)

        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    def as_python_source(self, tree_obj_name=None, tree_args=None, oids=False):
        """
        Returns string that will rebuild this tree in Python.
        """
        p = []

        if tree_obj_name is None:
            tree_obj_name = "tree_%s" % id(self)
        if self.label is not None:
            label = "'" + self.label + "'"
        else:
            label = "None"
        if oids:
            oid_str = ', oid="%s"' % self.oid
        else:
            oid_str = ""
        if tree_args is None:
            tree_args = ""
        else:
            tree_args = ", " + tree_args
        p.append("%s = dendropy.Tree(label=%s%s%s)" \
            % (tree_obj_name,
               label,
               oid_str,
               tree_args))
        if oids:
            p.append("%s.seed_node.oid = '%s'" % (tree_obj_name, self.seed_node.oid))

        taxon_obj_namer = lambda x: "tax_%s" % id(x)
        for taxon in self.taxon_set:
            tobj_name = taxon_obj_namer(taxon)
            if taxon.label is not None:
                label = "'" + taxon.label + "'"
            else:
                label = "None"
            if oids:
                oid_str = ', oid="%s"' % taxon.oid
            else:
                oid_str = ""
            p.append("%s = %s.taxon_set.require_taxon(label=%s%s)" \
                % (tobj_name,
                   tree_obj_name,
                   label,
                   oid_str))

        node_obj_namer = lambda x: "nd_%s" % id(x)
        for node in self.preorder_node_iter():
            for child in node.child_nodes():
                if node is self.seed_node:
                    nn = "%s.seed_node" % tree_obj_name
                else:
                    nn = node_obj_namer(node)
                if child.label is not None:
                    label = "'" + child.label + "'"
                else:
                    label = "None"
                if child.taxon is not None:
                    ct = taxon_obj_namer(child.taxon)
                else:
                    ct = "None"
                if oids:
                    oid_str = ', oid="%s"' % child.oid
                else:
                    oid_str = ""
                p.append("%s = %s.new_child(label=%s, taxon=%s, edge_length=%s%s)" %
                        (node_obj_namer(child),
                         nn,
                         label,
                         ct,
                         child.edge.length,
                         oid_str))
                if oids:
                    p.append('%s.edge.oid = "%s"' % (node_obj_namer(child), child.edge.oid))

        return "\n".join(p)

    def as_newick_string(self, **kwargs):
        """
        kwargs["reverse_translate"] can be function that takes a taxon and
        returns the label to appear in the tree.
        """
        return self.seed_node.as_newick_string(**kwargs)

    def print_newick(self, **kwargs):
        """
        Convenience method to newick string representation of this tree
        to the standard output stream.
        """
        import sys
        sys.stdout.write(self.as_newick_string(**kwargs))
        sys.stdout.write("\n")

    def as_ascii_plot(self, **kwargs):
        """
        Returns a string representation a graphic of this tree using ASCII
        characters.

        Keyword arguments:

            ``plot_metric``
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            ``show_internal_node_labels``
                Boolean: whether or not to write out internal node labels.
            - `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            ``leaf_spacing_factor``
                Positive integer: number of rows between each leaf.
            ``display_width``
                Force a particular display width, in terms of number of columns.

        """
        ap = AsciiTreePlot(**kwargs)
        return ap.compose(self)

    def write_ascii_plot(self, stream, **kwargs):
        """
        Writes an ASCII text graphic of this tree to `stream`.

        Keyword arguments:

            ``plot_metric``
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            ``show_internal_node_labels``
                Boolean: whether or not to write out internal node labels.
            - `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            ``leaf_spacing_factor``
                Positive integer: number of rows between each leaf.
            ``display_width``
                Force a particular display width, in terms of number of columns.

        """
        return stream.write(self.as_ascii_plot(**kwargs))

    def print_plot(self, **kwargs):
        """
        Writes an ASCII text graphic of this tree to standard output.

        Keyword arguments:

            ``plot_metric``
                A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            ``show_internal_node_labels``
                Boolean: whether or not to write out internal node labels.
            - `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            ``leaf_spacing_factor``
                Positive integer: number of rows between each leaf.
            ``display_width``
                Force a particular display width, in terms of number of columns.

        """
        import sys
        self.write_ascii_plot(sys.stdout, **kwargs)
        sys.stdout.write("\n")

    ###########################################################################
    ## Debugging/Testing

    def assign_node_labels_from_taxon_or_oid(self):
        for nd in self.postorder_node_iter():
            if nd.label is not None:
                continue
            if nd.taxon is not None:
                nd.label = nd.taxon.label
            else:
                nd.label = nd.oid

    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        if kwargs.get("splits"):
            if not kwargs.get("taxon_set"):
                kwargs["taxon_set"] = self.taxon_set
        self.seed_node.write_indented_form(out, **kwargs)

    def write_as_dot(self, out, **kwargs):
        """Writes the tree to `out` as a DOT formatted digraph"""
        if not kwargs.get("taxon_set"):
            kwargs["taxon_set"] = self.taxon_set
        out.write("digraph G {\n")

        nd_id_to_dot_nd = {}
        for n, nd in enumerate(self.preorder_node_iter()):
            label = format_node(nd, **kwargs)
            if nd is self.seed_node:
                label = "root %s" % label
            dot_nd = "n%d" % n
            out.write(' %s  [label="%s"];\n' % (dot_nd, label))
            nd_id_to_dot_nd[nd] = dot_nd
        for nd, dot_nd in nd_id_to_dot_nd.iteritems():
            try:
                e = nd.edge
                par_dot_nd = nd_id_to_dot_nd[e.tail_node]
            except:
                pass
            else:
                label = format_edge(e, **kwargs)
                s = ' %s -> %s [label="%s"];\n' % (par_dot_nd, dot_nd, label)
                out.write(s)
        out.write("}\n")

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
            `check_splits` if True specifies that the split_edge and split_bitmask attributes
                are checked.
        """
        check_splits = kwargs.get('check_splits', False)
        taxon_set = kwargs.get('taxon_set')
        if taxon_set is None:
            taxon_set = self.taxon_set
        if check_splits:
            taxa_mask = self.seed_node.edge.split_bitmask
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
                split_bitmask = curr_edge.split_bitmask
                assert((split_bitmask | taxa_mask) == taxa_mask)
            c = curr_node.child_nodes()
            if c:
                for child in c:
                    assert child.parent_node is curr_node
                    if check_splits:
                        cm |= child.edge.split_bitmask
            elif check_splits:
                assert(curr_node.taxon)
                cm = taxon_set.taxon_bitmask(curr_node.taxon)
            if check_splits:
                assert((cm & taxa_mask) == split_bitmask)
                assert self.split_edges[split_bitmask] == curr_edge
            curr_node, level = _preorder_list_manip(curr_node, siblings, ancestors)
        if check_splits:
            for s, e in self.split_edges.iteritems():
                assert(e in edges)
        return True

    def compose_newick(self):
        return self.as_newick_string(preserve_spaces=True)

##############################################################################
## Node

class Node(TaxonLinked):
    """
    A node on a tree, implementing only fundamental behaviour and
    properties.
    """

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

    ## INSTANCE METHODS #######################################################

    def __init__(self, **kwargs):
        TaxonLinked.__init__(self,
                             taxon=kwargs.get("taxon", None),
                             label=kwargs.get("label", None),
                             oid=kwargs.get("oid", None))
        self.age = None
        self._edge = None
        self._child_nodes = []
        self._parent_node = None
        self.edge = kwargs.get("edge", Edge(head_node=self))
        self._edge.head_node = self
        self.comments = []

    def __deepcopy__(self, memo):
        o = TaxonLinked.__deepcopy__(self, memo)
        for k, v in self.__dict__.iteritems():
            if not k in ['_child_nodes', '_taxon', "_oid"]:
                o.__dict__[k] = copy.deepcopy(v, memo)
        for c in self.child_nodes():
            o.add_child(copy.deepcopy(c, memo))
        memo[id(self._child_nodes)] = o._child_nodes
        memo[id(self._oid)] = o._oid
        return o

    ###########################################################################
    ## Iterators

    def preorder_iter(self, filter_fn=None):
        """
        Preorder traversal of self and its child_nodes.  Returns self
        and all descendants such that a node is returned before its
        child_nodes (and their child_nodes). Filtered by filter_fn: node is
        only returned if no filter_fn is given or if filter_fn returns
        True.
        """
        stack = [self]
        while stack:
            node = stack.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            child_nodes.extend(stack)
            stack = child_nodes

    def postorder_iter(self, filter_fn=None):
        """
        Postorder traversal of the self and its child_nodes.  Returns self
        and all descendants such that a node's child_nodes (and their
        child_nodes) are visited before node.  Filtered by filter_fn:
        node is only returned if no filter_fn is given or if filter_fn
        returns True.
        """
        stack = [(self, False)]
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

    def leaf_iter(self, filter_fn=None):
        """
        Returns an iterator over the leaf_nodes that are descendants of self
        (with leaves returned in same order as a post-order traversal of the tree).
        """
        if filter_fn:
            filter_fn = lambda x: x.is_leaf() and filter_fn(x) or None
        else:
            filter_fn = lambda x: x.is_leaf() and x or None
        for node in self.postorder_iter(filter_fn):
            yield node

    def level_order_iter(self, filter_fn=None):
        """
        Level-order traversal of self and its child_nodes. Filtered
        by filter_fn: node is only returned if no filter_fn is given
        or if filter_fn returns True
        """
        if filter_fn is None or filter_fn(self):
            yield self
        remaining = self.child_nodes()
        while len(remaining) > 0:
            node = remaining.pop(0)
            if filter_fn is None or filter_fn(node):
                yield node
            child_nodes = node.child_nodes()
            remaining.extend(child_nodes)

    def ancestor_iter(self, filter_fn=None, inclusive=True):
        """
        Iterates over all ancestors of self. If `inclusive` is True, self
        is returned as the first item of the sequence.
        """
        if inclusive:
            yield self
        node = self
        while node is not None:
            node = node.parent_node
            if node is not None \
                   and (filter_fn is None or filter_fn(node)):
                yield node

    def age_order_iter(self, include_leaves=True, filter_fn=None, descending=False):
        """
        Iterates over nodes in order of age. If `include_leaves` is False, will
        skip leaves (default is not to skip leaves). If `descending` is True,
        will go from oldest nodes to youngest (default is asecending: youngest
        nodes to oldest).
        """
        if not descending:
            leaves = [nd for nd in self.leaf_iter()]
            queued_pairs = []
            in_queue = set()
            for leaf in leaves:
                age_nd_tuple = (leaf.age, leaf)
                queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
                in_queue.add(leaf)
            while queued_pairs:
                next_el = queued_pairs.pop(0)
                age, nd = next_el
                in_queue.remove(nd)
                p = nd.parent_node
                if p and p not in in_queue:
                    age_nd_tuple = (p.age, p)
                    queued_pairs.insert(bisect.bisect(queued_pairs, age_nd_tuple), age_nd_tuple)
                    in_queue.add(p)
                if include_leaves or nd.is_internal():
                    yield nd
        else:
            nds = [(nd.age, nd) for nd in self.preorder_iter()]
            nds.sort(reverse=True)
            for nd in nds:
                if include_leaves or nd[1].is_internal():
                    yield nd[1]

    def postorder_internal_node_iter(self, filter_fn=None):
        """
        Iterates over all internal nodes in post-order.
        """
        if filter_fn:
            filter_fn = lambda x: (not x.is_leaf() and filter_fn(x)) or None
        else:
            filter_fn = lambda x: (x and not x.is_leaf()) or None
        for node in self.postorder_iter(filter_fn):
            yield node

    def preorder_internal_node_iter(self, filter_fn=None):
        """
        Iterates over all internal nodes in pre-order.
        """
        if filter_fn:
            filter_fn = lambda x: (not x.is_leaf() and filter_fn(x)) or None
        else:
            filter_fn = lambda x: (x and not x.is_leaf()) or None
        for node in self.preorder_iter(filter_fn):
            yield node

    ###########################################################################
    ## (Attribute) Accessors and Mutators

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

    ###########################################################################
    ## Metrics

    def level(self):
        "Number of nodes between self and root."
        if self.parent_node:
            return self.parent_node.level() + 1
        else:
            return 0

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

    def distance_from_tip(self):
        """
        Sum of edge lengths from tip to node. If tree is not ultrametric
        (i.e., descendent edges have different lengths), then count the
        maximum of edge lengths. Note that the 'calc_node_ages()' method
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

    ###########################################################################
    ## Structural Infornation and Manipulation

    def leaf_nodes(self):
        """
        Returns list of all leaf_nodes descended from this node (or just
        list with self as the only member if self is a leaf).
        """
        return [node for node in \
                self.postorder_iter(lambda x: bool(len(x.child_nodes())==0))]

    def child_nodes(self):
        "Returns the a shallow-copy list of all child nodes."
        return list(self._child_nodes)

    def set_child_nodes(self, child_nodes):
        """
        Sets the child_nodes for this node.
        Side effects:

            - sets the parent of each child node to this node
            - sets the tail node of each child to self
        """
        self._child_nodes = list(child_nodes)
        for nidx in range(len(self._child_nodes)):
            self._child_nodes[nidx].parent_node = self
            self._child_nodes[nidx].edge.tail_node = self

    def set_children(self, child_nodes):
        """Legacy support: delegates to `set_child_nodes()`"""
        return self.set_child_nodes(child_nodes)

    def _get_parent_node(self):
        """Returns the parent node of this node."""
        return self._parent_node

    def _set_parent_node(self, parent):
        """Sets the parent node of this node."""
        self._parent_node = parent
        self.edge.tail_node = parent

    parent_node = property(_get_parent_node, _set_parent_node)

    def incident_edges(self):
        """Return parent and child edges."""
        e = [c.edge for c in self._child_nodes]
        e.append(self.edge)
        return e

    def get_incident_edges(self):
        """Legacy synonym for 'incident_edges()'"""
        return self.incident_edges()

    def adjacent_nodes(self):
        """Return parent and child nodes."""
        n = [c for c in self._child_nodes]
        if self.parent_node:
            n.append(self.parent_node)
        return n

    def get_adjacent_nodes(self):
        """Legacy synonym for 'get_incident_edges()'"""
        return self.adjacent_nodes()

    def sister_nodes(self):
        """Return all other children of parent, excluding self."""
        p = self.parent_node
        if not p:
            return []
        sisters = [nd for nd in p.child_nodes() if nd is not self]
        return sisters

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
        """This function should be used to "undo" the effects of
        Node.reversible_remove_child NOTE: the behavior is only
        guaranteed if the tree has not been modified between the
        remove_child and reinsert_nodes calls! (or the tree has been
        restored such that the node/edge identities are identical to the
        state before the remove_child call.

        The order of info in each tuple is:

            0 - node removed
            1 - parent of node removed
            2 - pos in parent matrix
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

    def collapse_clade(self):
        """Collapses all internal edges that are descendants of self."""
        if self.is_leaf():
            return
        leaves = [i for i in self.leaf_iter()]
        self.set_child_nodes(leaves)

    ###########################################################################
    ## Representation

    def description(self, depth=1, indent=0, itemize="", output=None, taxon_set=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sNode object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Edge]' % leader1)
            if self.edge is not None:
                edge_desc = self.edge.description(0)
            else:
                edge_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, edge_desc))

            output_strio.write('\n%s[Taxon]' % leader1)
            if self.taxon is not None:
                taxon_desc = self.taxon.description(0)
            else:
                taxon_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, taxon_desc))

            output_strio.write('\n%s[Parent]' % leader1)
            if self.parent_node is not None:
                parent_node_desc = self.parent_node.description(0)
            else:
                parent_node_desc = 'None'
            output_strio.write('\n%s%s' % (leader2, parent_node_desc))
            if hasattr(self.edge, 'split_bitmask'):
                output_strio.write('\n%s[Clade Mask]' % leader1)
                if taxon_set is None:
                    output_strio.write('\n%s%s' % (leader2, self.edge.split_bitmask))
                else:
                    output_strio.write('\n%s%s' % (leader2, taxon_set.split_bitmask_string(self.edge.split_bitmask)))
            output_strio.write('\n%s[Children]' % leader1)
            if len(self._child_nodes) == 0:
                output_strio.write('\n%sNone' % leader2)
            else:
                for i, cnd in enumerate(self._child_nodes):
                    output_strio.write('\n%s[%d] %s' % (leader2, i, cnd.description(0)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ###########################################################################
    ## For debugging we build-in a full-fledged NEWICK composition independent
    ## of the nexus/newick family of modules. Client code should prefer to
    ## use Newick/Nexus readers/writers, or Tree.write(), TreeList.write(),
    ## DataSet.write() etc.

    def as_newick_string(self, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.
        For production purposes, use the the full-fledged 'as_string()'
        method of the object.
        """
        out = StringIO()
        self.write_newick(out, **kwargs)
        return out.getvalue()

    def write_newick(self, out, **kwargs):
        """
        This returns the Node as a NEWICK statement according to the given
        formatting rules. This should be used for debugging purposes only.  For
        production purposes, use the the full-fledged 'write_to_stream()'
        method of the object.
        """
        edge_lengths = not kwargs.get('suppress_edge_lengths', False)
        edge_lengths = kwargs.get('edge_lengths', edge_lengths)
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
        `suppress_internal_labels` is a Boolean, and defaults to False.
        """
        is_leaf = (len(self._child_nodes) == 0)
        if not is_leaf:
            if kwargs.get("newick", False):
                return self.as_newick_string()
            if kwargs.get("suppress_internal_labels", False) \
                    or not kwargs.get("include_internal_labels", True):
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
        preserve_spaces = kwargs.get("preserve_spaces", False)
        raw_labels = kwargs.get("raw_labels", False)
        if raw_labels:
            return tag
        return textutils.escape_nexus_token(tag, preserve_spaces=preserve_spaces)

    ###########################################################################
    ## alternate representation of tree structure for debugging

    def get_indented_form(self, **kwargs):
        out = StringIO()
        self.write_indented_form(out, **kwargs)
        return out.getvalue()

    def write_indented_form(self, out, **kwargs):
        indentation = kwargs.get("indentation", "    ")
        split_bitmasks = kwargs.get("splits", True)
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
            cm = "%s " % format_split(self.edge.split_bitmask, **kwargs)
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
        __init__ creates an edge from tail_node to head_node.  Modified from
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

    ###########################################################################
    ## Structural

    def collapse(self):
        """
        Inserts all children of the head_node of self as children of the
        tail_node of self in the same place in the child_node list that head_node
        had occupied. The edge length and head_node will no longer be
        part of the tree.
        """
        to_del = self.head_node
        p = self.tail_node
        if not p:
            return
        children = to_del.child_nodes()
        if not children:
            raise ValueError('collapse_self called with a terminal.')
        pos = p.child_nodes().index(to_del)
        p.remove_child(to_del)
        for child in children:
            p.add_child(child, pos=pos)
            pos += 1

    ###########################################################################
    ## Representation

    def description(self, depth=1, indent=0, itemize="", output=None, taxon_set=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return
        output_strio = StringIO()
        if self.label is None:
            label = " (%s, Length=%s)" % (self.oid, str(self.length))
        else:
            label = " (%s: '%s', Length=%s)" % (self.oid, self.label, str(self.length))
        output_strio.write('%s%sEdge object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            leader1 = ' ' * (indent + 4)
            leader2 = ' ' * (indent + 8)
            output_strio.write('\n%s[Length]' % leader1)
            if self.length is not None:
                length = self.length
            else:
                length = "None"
            output_strio.write('\n%s%s' % (leader2, length))
            output_strio.write('\n%s[Tail Node]' % leader1)
            if self.tail_node is not None:
                tn = self.tail_node.description(0)
            else:
                tn = "None"
            output_strio.write('\n%s%s' % (leader2, tn))
            output_strio.write('\n%s[Head Node]' % leader1)
            if self.head_node is not None:
                hn = self.head_node.description(0)
            else:
                hn = "None"
            output_strio.write('\n%s%s' % (leader2, hn))
            if hasattr(self, 'split_bitmask'):
                output_strio.write('\n%s[Clade Mask]' % leader1)
                if taxon_set is None:
                    output_strio.write('\n%s%s' % (leader2, self.split_bitmask))
                else:
                    output_strio.write('\n%s%s' % (leader2, taxon_set.split_bitmask_string(self.split_bitmask)))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

###############################################################################
## AsciiTreePlot

class AsciiTreePlot(object):

    class NullEdgeLengthError(ValueError):
        def __init__(self, *args, **kwargs):
            ValueError.__init__(self, *args, **kwargs)

    def __init__(self, **kwargs):
        """
        __init__ takes the following kwargs:

            - `plot_metric` A string which specifies how branches should be scaled, one of:
                'age' (distance from tips), 'depth' (distance from root),
                'level' (number of branches from root) or 'length' (edge
                length/weights).
            - `show_internal_node_labels`
                Boolean: whether or not to write out internal node labels.
            - `show_internal_node_ids`
                Boolean: whether or not to write out internal node id's.
            - `leaf_spacing_factor`
                Positive integer: number of rows between each leaf.
            - `display_width`
                Force a particular display width, in terms of number of columns.

        """
        self.plot_metric = kwargs.get('plot_metric', 'depth')
        self.show_internal_node_labels = kwargs.get('show_internal_node_labels', False)
        self.show_internal_node_ids = kwargs.get('show_internal_node_ids', False)
        self.leaf_spacing_factor = kwargs.get('leaf_spacing_factor', 2)
#        self.null_edge_length = kwargs.get('null_edge_length', 0)
        self.display_width = kwargs.get('display_width', None)
        self.reset()

    def reset(self):
        self.grid = []
        self.node_row = {}
        self.node_col = {}
        self.node_offset = {}
        self.current_leaf_row = 0

    def _calc_node_offsets(self, tree):
        if self.plot_metric == 'age' or self.plot_metric == 'depth':

            ## for verification ...
#            tree.calc_node_ages(check_prec=False)
#            for nd in tree.postorder_node_iter():
#                self.node_offset[nd] = nd.age
#            flipped_origin = max(self.node_offset.values())
#            for nd in self.node_offset:
#                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
#            return

            for nd in tree.postorder_node_iter():
                cnds = nd.child_nodes()
                if self.plot_metric == 'depth': # 'number of branchings from tip'
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        depths = [self.node_offset[v] for v in cnds]
                        curr_node_offset = max(depths) + 1
#                        print curr_node_offset, [self.node_offset[v] for v in cnds]
                elif self.plot_metric == 'age': # 'sum of edge weights from tip'
                    # note: no enforcement of ultrametricity!
                    if len(cnds) == 0:
                        curr_node_offset = 0.0
                    else:
                        if cnds[0].edge.length is not None:
                            curr_node_offset = self.node_offset[cnds[0]] + cnds[0].edge.length
#                        if len(elens) == 0:
#                            curr_node_offset = self.node_offset[cnds[0]]
#                        else:
#                            curr_node_offset = max(elens) + self.node_offset[cnds[0]]
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                self.node_offset[nd] = curr_node_offset
            flipped_origin = max(self.node_offset.values())
            for nd in self.node_offset:
                self.node_offset[nd] = flipped_origin - self.node_offset[nd]
        else:
            for nd in tree.preorder_node_iter():
                if self.plot_metric == 'level': # 'number of branchings from root'
                    curr_edge_len = 1
                elif self.plot_metric == 'length': # 'sum of edge weights from root'
                    if nd.edge.length is not None:
                        curr_edge_len = nd.edge.length
                    else:
                        curr_edge_len = 0
                else:
                    raise ValueError("Unrecognized plot metric '%s' (must be one of: 'age', 'depth', 'level', or 'length')" % self.plot_metric)
                if nd.parent_node is None:
                    self.node_offset[nd] = curr_edge_len
                else:
                    self.node_offset[nd] =  curr_edge_len + self.node_offset[nd.parent_node]
#        print "\n".join([str(k) for k in self.node_offset.values()])

    def draw(self, tree, dest):
        dest.write(self.compose(tree))

    def get_label_for_node(self, nd):
        if nd.taxon and nd.taxon.label:
            return nd.taxon.label
        # @TODO: we should have a separate setting for labeling nodes with an
        # id, but thus far when I want to see this, I want
        # internal_nodes_labels too...
        label = []
        if self.show_internal_node_labels and nd.label:
            label.append(nd.label)
        if self.show_internal_node_ids:
            label.append("@")
            label.append(str(id(nd)))
        if not label:
            return "@"
        return "".join(label)

    def compose(self, tree):
        self.reset()
        if self.display_width is None:
            display_width = termutils.terminal_width() - 1
        else:
            display_width = self.display_width
        max_label_len = max([len(self.get_label_for_node(i)) for i in tree.leaf_iter()])
        if max_label_len <= 0:
            max_label_len = 0
        #effective_display_width = display_width - max_label_len - len(tree.internal_nodes) - 1
        effective_display_width = display_width - max_label_len - 1
        self._calc_node_offsets(tree)
        widths = [self.node_offset[i] for i in tree.leaf_iter() if self.node_offset[i] is not None]
        max_width = float(max(widths))
        if max_width == 0:
            raise AsciiTreePlot.NullEdgeLengthError("Tree cannot be plotted under metric '%s' due to zero or null edge lengths: '%s'" % (self.plot_metric, tree.as_newick_string()))
        edge_scale_factor = float(effective_display_width) / max_width
        self.calc_plot(tree.seed_node,
                       edge_scale_factor=edge_scale_factor)
        for i in range(len(tree.leaf_nodes())*self.leaf_spacing_factor + 1):
            self.grid.append([' ' for i in range(0, display_width)])
        self.draw_node(tree.seed_node)
        display = '\n'.join([''.join(i) for i in self.grid])
        return display

    def calc_plot(self, node, edge_scale_factor):
        """
        First pass through tree, post-order traversal to calculate
        coordinates of each node.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for n in child_nodes:
                self.calc_plot(n, edge_scale_factor)
            ys = [self.node_row[n] for n in child_nodes]
            self.node_row[node] = int(float((max(ys)-min(ys)) / 2) + min(ys))
        else:
            self.node_row[node] = self.current_leaf_row
            self.current_leaf_row = self.current_leaf_row + self.leaf_spacing_factor
        if node.edge.length is None:
            self.node_col[node] = 1
        else:
            self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)
        self.node_col[node] = int(float(self.node_offset[node]) * edge_scale_factor)

    def draw_label(self, label, row, start_col):
        if label:
            for i in range(len(label)):
                if start_col + i < len(self.grid[row]):
                    self.grid[row][start_col+i] = label[i]

    def draw_node(self, node):
        """
        Second pass through tree, plotting nodes onto given self.grid.
        """
        child_nodes = node.child_nodes()
        if child_nodes:
            for i, child_node in enumerate(child_nodes):
                start_row = min([self.node_row[node], self.node_row[child_node]])
                end_row = max([self.node_row[node], self.node_row[child_node]])
                if i == 0:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '/'
                    start_row = start_row+1
                    edge_row = self.node_row[child_node]
                elif i == len(child_nodes)-1:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '\\'
                    edge_row = self.node_row[child_node]
                else:
                    self.grid[self.node_row[child_node]][self.node_col[node]] = '+'
                    edge_row = self.node_row[child_node]
                self.draw_node(child_node)
                for x in range(self.node_col[node]+1, self.node_col[child_node]):
                    self.grid[edge_row][x] = '-'
                for y in range(start_row, end_row):
                    self.grid[y][self.node_col[node]] = '|'
            label = []
            if self.show_internal_node_labels or self.show_internal_node_ids:
                label = self.get_label_for_node(node)
                self.draw_internal_text(label, self.node_row[node], self.node_col[node])
            else:
                self.grid[self.node_row[node]][self.node_col[node]]='+'
        else:
            label = self.get_label_for_node(node)
            self.draw_label(label, self.node_row[node], self.node_col[node]+1)

    def draw_internal_text(self, label, r, c):
        row = self.grid[r]
        try:
            for n, letter in enumerate(label):
                row[c + n] = letter
        except:
            pass

###############################################################################
## NodeRelationship

class NodeRelationship(object):

    def from_tree(tree):
        return [NodeRelationship(node) for node in tree]
    from_tree = staticmethod(from_tree)

    def from_node(node):
        ndr = NodeRelationship(None, None, None, None)
        if node.parent_node is not None:
            ndr.parent_label = node.parent_node.label
        else:
            ndr.parent_label = 'None'
        ndr.child_labels = [cnd.label for cnd in node.child_nodes()]
        ndr.edge_length = node.edge.length
        if node.taxon is not None:
            ndr.taxon_label = node.taxon.label
        else:
            ndr.taxon_label = 'None'
        return ndr
    from_node = staticmethod(from_node)

    def __init__(self, parent_label, child_labels, edge_length, taxon_label):
        self.parent_label = parent_label
        self.child_labels = child_labels
        self.edge_length = edge_length
        self.taxon_label = taxon_label

    def test_node(self, testcase, node):
        if self.parent_label is not None:
            testcase.assertTrue(node.parent_node is not None)
            testcase.assertEqual(self.parent_label, node.parent_node.label)
        else:
            testcase.assertTrue(node.parent_node is None or node.parent_node.label is None)
        testcase.assertEqual(self.child_labels, [cnd.label for cnd in node.child_nodes()])
        if self.edge_length is not None:
            testcase.assertTrue(node.edge.length is not None)
            testcase.assertAlmostEqual(self.edge_length, node.edge.length)
        else:
            testcase.assertTrue(node.edge.length is None)
        if self.taxon_label is not None:
            testcase.assertTrue(node.taxon is not None)
            testcase.assertEqual(self.taxon_label, node.taxon.label)
        else:
            testcase.assertTrue(node.taxon is None or node.taxon.label is None)

###############################################################################
## Helper Functions

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
    nf = kwargs.get('node_formatter', None)
    if nf:
        return nf(nd)
    if nd.is_leaf():
        t = nd.taxon
        if t:
            label = t.label
        else:
            label = "anonymous leaf"
    else:
        label = "* %s" % str(nd.oid)
    return label

def format_edge(e, **kwargs):
    ef = kwargs.get('edge_formatter', None)
    if ef:
        return ef(e)
    return str(e)

def format_split(split, width=None, **kwargs):
    from dendropy.treesplit import split_as_string
    if width is None:
        width = len(kwargs.get("taxon_set"))
    s = split_as_string(split, width, symbol1=kwargs.get("off_symbol"), symbol2=kwargs.get("on_symbol"))
    return s

def convert_node_to_root_polytomy(nd):
    """If `nd` has two children and at least on of them is an internal node,
    then it will be converted to an out-degree three node (with the edge length
    added as needed).

    Returns a tuple of child nodes that were detached (or() if the tree was not
    modified). This can be useful for removing the deleted node from the split_edges
    dictionary.
    """
    nd_children = nd.child_nodes()
    if len(nd_children) > 2:
        return ()
    try:
        left_child = nd_children[0]
    except:
        return ()
    if not left_child:
        return ()
    if len(nd_children) == 1:
        right_child = None
        dest_edge_head = nd
    else:
        right_child = nd_children[1]
        dest_edge_head = right_child
    curr_add = None
    if right_child and right_child.is_internal():
        try:
            left_child.edge.length += right_child.edge.length
        except:
            pass
        nd.remove_child(right_child)
        grand_kids = right_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = right_child
    elif left_child.is_internal():
        try:
            dest_edge_head.edge.length += left_child.edge.length
        except:
            pass
        nd.remove_child(left_child)
        grand_kids = left_child.child_nodes()
        for gc in grand_kids:
            nd.add_child(gc)
        curr_add = left_child
    if curr_add:
        ndl = [curr_add]
        t = convert_node_to_root_polytomy(nd)
        ndl.extend(t)
        return tuple(ndl)

    return ()
