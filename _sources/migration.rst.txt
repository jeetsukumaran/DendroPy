#######################################
DendroPy 4 Changes and Migration Primer
#######################################

Introduction
============

* Updated for full (and exclusive) Python 3.x compatibility.

* Faster, better, stronger! Core serialization/deserialization infrastructure
  rewritten from the ground up, with *many* optimizations for speed and
  reliability.

Python Version Compatibility
============================

* Compatibility: Python 3 is fully supported. The only version of Python 2
  supported is Python 2.7.

    * Python 2: Python 2.7

    * Python 3: Python 3.1, 3.2, 3.3, 3.4

Library-Wide Changes
====================

Public Module Reorganization
----------------------------

A number of modules have been renamed, moved, or split into multiple modules.
Calls to the old module should continue to work, albeit with warnings exhorting
that you update to the latest configuration.

    * ``dendropy.treecalc`` has been split into three submodules depending on
      whether the statistic or value being calculated is on a single tree, a
      single tree and a dataset, or two trees:

        *   :mod:`dendropy.calculate.treemeasure`
            For calculation of statistics, metrics, and values on a single tree.
        *   :mod:`dendropy.calculate.treecompare`
            For calculation of statistics, metrics, and values of two trees
            (e.g., Robinson-Fould's distances).
        *   :mod:`dendropy.calculate.treescore`
            For calculation of statistics, metrics, and values of a tree with
            reference to a dataset under some criterion.
    * The functionality provided ``dendropy.treesplit`` has been largely subsumed by the new |Bipartition| class.
    * The functionality provided by ``dendropy.treesum`` has been largely subsumed by the new |TreeArray| class, a high-performance class for efficiently managing and carrying out operations on large collections of large trees.
    * ``dendropy.reconcile`` has been moved to :mod:`dendropy.model.reconcile`.
    * ``dendropy.coalescent`` has been moved to :mod:`dendropy.model.coalescent`.
    * ``dendropy.popgenstat`` has been moved to :mod:`dendropy.calculate.popgenstat`.
    * ``dendropy.treesim`` has been moved to :mod:`dendropy.simulate.treesim`.
    * ``dendropy.popgensim`` has been moved to :mod:`dendropy.simulate.popgensim`.


Behind-the-Scenes Module Reorganization
---------------------------------------

* In constrast to the above, the following changes *should* be opaque to most
  normal usage and client code. Most of the names (classes/methods/variables)
  in these modules were imported into the '``dendropy``' namespace, and this is
  how all public code should be accessing them, *or* they were never exposed
  (or meant to be exposed) for public usage in the first place. A list of
  module changes:

        +----------------------------------+-----------------------------------------------+
        | DendroPy 3                       | DendroPy 4                                    |
        +==================================+===============================================+
        | :mod:`dendropy.dataobject.base`  | :mod:`dendropy.datamodel.basemodel`           |
        +----------------------------------+-----------------------------------------------+
        | :mod:`dendropy.dataobject.taxon` | :mod:`dendropy.datamodel.taxonmodel`          |
        +----------------------------------+-----------------------------------------------+
        | :mod:`dendropy.dataobject.tree`  | :mod:`dendropy.datamodel.treemodel`           |
        |                                  | :mod:`dendropy.datamodel.treecollectionmodel` |
        +----------------------------------+-----------------------------------------------+
        | :mod:`dendropy.dataobject.char`  | :mod:`dendropy.datamodel.charstatemodel`,     |
        |                                  | :mod:`dendropy.datamodel.charmatrixmodel`     |
        +----------------------------------+-----------------------------------------------+


Unique Object Identifier ("``oid``") Attributes Removed
-------------------------------------------------------

* The entire ``oid`` system ("object identifier"), i.e., the unique id assigned
  to every data object, has been removed. This was an implementation artifact
  from NEXML parsing that greatly slowed down a number of operations without
  any benefit or utility for most normal operations.

:class:`TaxonSet` is now :class:`TaxonNamespace`
================================================

* The ``dendropy.TaxonSet`` class has been renamed |TaxonNamespace|,
  (and the corresponding ``taxon_set`` attribute of phylogenetic data objects
  that reference a taxonomic context has been renamed ``taxon_namespace``).

* The |TaxonNamespace| class replaces the :class:`TaxonSet` class as the
  manager for the :class:`Taxon` objects.

* The API is largely similar with the following differences:

    * Calls to the
      :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.__getitem__` and
      :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.__delitem__` methods
      (e.g. ``TaxonNamespace[x]``) now only accept integer values as arguments
      (representing indexes into the list of :class:`Taxon` objects in the
      internal array).

    * :meth:`TaxonSet.has_taxon()` and :meth:`TaxonSet.has_taxa()` have been
        replaced by :meth:`TaxonNamespace.has_taxon_label()` and
        :meth:`TaxonNamespace.has_taxa_labels()` respectively.

    * Various new methods for accessing and managing the collection of
        :class:`Taxon` objects (e.g., :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.findall`, :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.remove_taxon`, :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.remove_taxon_label`,
        :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.discard_taxon_label`, :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespace.__delitem__`, etc.)

    * Numerous look-up methods took '``case_insensitive``' as an argument that
      determined whether the look-up was case sensitive or not (when
      retrieving, for example, a :class:`Taxon` object corresponding to a
      particular label), which, if not specified, default to ``False``, i.e. a
      non-caseless or a case-sensitive matching criteria. In all cases, this
      has been changed to to '``case_sensitive``' with a default of ``True``. That
      is, searches by default are still case-sensitive, but now you will have
      to specify '``case_sensitive=False``' instead of '``case_insensitive=True``'
      to perform a case-insensitive search. This change was for consistency
      with the rest of the library.

* In most cases, a simple global search-and-replace of "TaxonSet" with
  "TaxonNamespace" and "``taxon_set``" with "``taxon_namespace``" should be
  sufficient to bring existing code into line with DendroPy 4.

* For legacy support, a class called :class:`TaxonSet` exists. This derives with no
  modifications from :class:`TaxonNamespace`. Instantiating objects of this class
  will result in warnings being emitted. As long as usage of :class:`TaxonSet` does
  conforms to the above API change notes, old or legacy code should continue
  to work unchanged (albeit, with some warning noise). This support is
  temporary and will be removed in upcoming releases: code should update to
  using :class:`TaxonNamespace` as soon as expedient.

* For legacy support, "``taxon_set``" continues to be accepted and processed as
  an attribute name and keyword argument synonymous with "``taxon_namespace``".
  Usage of this will result in warnings being emitted, but code should
  continue to function as expected. This support is temporary and will be
  removed in upcoming releases: code should update to using
  "``taxon_namespace``" as soon as expedient.

The :class:`Node` Class
=======================

* Constructor now only accepts keyword arguments (and ``oid`` is *not* one of them!).

* :meth:`~dendropy.datamodel.treemodel.Node.add_child()` no longer accepts ``pos`` as an argument to indicate
  position in which a child should be inserted. Use :meth:`~dendropy.datamodel.treemodel.Node.insert_child()`
  which takes a position specified by ``index`` and a node specified by ``node``
  for this functionality instead.

The :class:`Edge` Class
=======================

* Constructor now only accepts keyword arguments (and ``oid`` is *not* one of them!).

* Because ``tail_node`` is no longer an independent attribute but a dynamic
  property, bound to :attr:`Node._parent_node` attribute of the ``head_node``
  (see below), the :class:`Edge` constructor does *not* accept ``tail_node`` as
  an argument.

* The ``tail_node`` of an :class:`Edge` object is now a dynamic property,
  referencing the :attr:`Node._parent_node` attribute of the
  :attr:`Edge._head_node` of the :class:`Edge` object. So, now updating
  :attr:`Edge._tail_node` of an :class:`Edge` object will set the
  :attr:`Node._parent_node` of its :attr:`Edge._head_node` to the new value,
  and vice versa.  This avoids the need for independent book-keeping logic to
  ensure that :attr:`Node._parent_node` and :attr:`Edge._tail_node` are always
  synchronized to reference the same :class:`Node` object and all the potential
  errors this might cause.

The :class:`Tree` Class
=======================

* Constructor no longer supports they ``stream`` keyword argument to construct
  the new :class:`~dendropy.datamodel.treemodel.Tree` object from a data source. Use the factory class
  method: :meth:`~dendropy.datamodel.treemodel.Tree.get_from_stream()` instead.

* :meth:`~dendropy.datamodel.treemodel.Tree.nodes()` : sorting option removed; use :func:`~dendropy.datamodel.treemodel.sorted(tree.nodes())` instead.

* :meth:`~dendropy.datamodel.treemodel.Tree.node_set()` : removed; use :func:`~dendropy.datamodel.treemodel.set(tree.nodes())` instead.

* :meth:`~dendropy.datamodel.treemodel.Tree.edge_set()` : removed; use :func:`~dendropy.datamodel.treemodel.set(tree.edges())` instead.

* For consistency with :meth:`~dendropy.datamodel.treemodel.Tree.preorder_node_iter()`,
  :meth:`~dendropy.datamodel.treemodel.Tree.postorder_node_iter()`, a number of iteration methods have been renamed.

    +----------------------------------+-------------------------------------------------------------------+
    | DendroPy 3                       | DendroPy 4                                                        |
    +==================================+===================================================================+
    | ``Tree.level_order_node_iter()`` | :meth:`~dendropy.datamodel.treemodel.Tree.levelorder_node_iter()` |
    +----------------------------------+-------------------------------------------------------------------+
    | ``Tree.level_order_edge_iter()`` | :meth:`~dendropy.datamodel.treemodel.Tree.levelorder_edge_iter()` |
    +----------------------------------+-------------------------------------------------------------------+
    | ``Node.level_order_iter()``      | :meth:`~dendropy.datamodel.treemodel.Node.levelorder_iter()`      |
    +----------------------------------+-------------------------------------------------------------------+
    | ``Tree.age_order_node_iter()``   | :meth:`~dendropy.datamodel.treemodel.Tree.ageorder_node_iter()`   |
    +----------------------------------+-------------------------------------------------------------------+
    | ``Node.age_order_iter()``        | :meth:`~dendropy.datamodel.treemodel.Node.ageorder_iter()`        |
    +----------------------------------+-------------------------------------------------------------------+
    | ``Tree.leaf_iter()``             | :meth:`~dendropy.datamodel.treemodel.Tree.leaf_node_iter()`       |
    +----------------------------------+-------------------------------------------------------------------+

  The old names are still supported for now (with warnings being emitted),
  but new code should start using the newer names.  In additon, support for
  in-order or infix tree traversal has been added:
  :meth:`~dendropy.datamodel.treemodel.Tree.inorder_node_iter`, :meth:`~dendropy.datamodel.treemodel.Tree.inorder_edge_iter()`.

* Instead of ``tree_source_iter`` and ``multi_tree_source_iter``, use
  :meth:`~dendropy.datamodel.treemodel.Tree.yield_from_files`

NEWICK-format Reading
=====================

* The ``suppress_external_taxon_labels`` and ``suppress_external_node_labels`` keyword
  arguments have been replaced by ``suppress_leaf_taxon_labels`` and
  ``suppress_leaf_node_labels``, respectively. This is for consistency with the
  rest of the library (including writing in NEWICK-format), which uses the term
  "leaf" rather than "external".

* The various boolean rooting directive switches (``as_rooted``,
  ``default_as_rooted``, etc.) have been replaced by a single argument:
  ``rooting``. This can take on one of the following (string) values:

    * rooting="default-unrooted"
        Interpret trees following rooting token ("``[&R]``" for rooted,
        "``[&U]``" for unrooted) if present; otherwise, intrepret trees as
        unrooted.
    * rooting"default-rooted"
        Interpret trees following rooting token ("``[&R]``" for rooted,
        "``[&U]``" for unrooted) if present; otherwise, intrepret trees as
        rooted.
    * rooting="force-unrooted"
        Unconditionally interpret all trees as unrooted.
    * rooting="force-rooted"
        Unconditionally interpret all trees as rooted.

  The value of the "``rooting``" argument defaults to "default-unrooted", i.e.,
  all trees are assumed to be unrooted unless a rooting token is present that
  explicitly specifies the rooting state.

NEWICK-format Writing
=====================

* Previously, if ``annotations_as_nhx`` was ``True``, metadata annotations would
  be written out even if ``suppress_annotations`` was ``True``. Now,
  ``suppress_annotations`` must be ``True`` for annotations to be written out,
  even if ``annotations_as_nhx`` is ``True``.

The :class:`DataSet` Class
==========================

* Constructor no longer supports they ``stream`` keyword argument to construct
  the new :class:`DataSet` object from a data source. Use the factory class
  method: :meth:`DataSet.get_from_stream()` instead.

* Constructor only accepts one unnamed (positional) argument: either a
  :class:`DataSet` instance to be cloned, or an iterable of
  :class:`TaxonNamespace`, :class:`TreeList`, or
  :class:`CharacterMatrix`-derived instances to be composed (added) into the
  new :class:`DataSet` instance.

* :class:`TaxonNamespace` no longer managed.



