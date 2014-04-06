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

Unique Object Identifier ("`oid`") Attributes Removed
=====================================================

* The entire `oid` system ("object identifier"), i.e., the unique id assigned
  to every data object, has been removed. This was an implementation artifact
  from NEXML parsing.


:class:`TaxonSet` is now :class:`TaxonNamespace`
================================================

* The `dendropy.TaxonSet` class has been renamed `dendropy.TaxonNamespace`
  (and the corresponding `taxon_set` attribute of phylogenetic data objects
  that reference a taxonomic context has been renamed `taxon_namespace`).

* The :class:`TaxonNamespace` class replaces the :class:`TaxonSet` class as the
  manager for the :class:`Taxon` objects.

* The API is largely similar with the following differences:

  * Calls to the `__getitem__()` and `__delitem()__` methods (e.g.
    'TaxonNamespace[x]') now only accept integer values as arguments
    (representing indexes into the list of :class:`Taxon` objects in the
    internal array).

  * :meth:`TaxonSet.has_taxon()` and :meth:`TaxonSet.has_taxa()` have been
    replaced by :meth:`TaxonNamespace.has_taxon_label()` and
    :meth:`TaxonNamespace.has_taxa_labels()` respectively.

  * Various new methods for accessing and managing the collection of
    :class:`Taxon` objects (e.g., `findall`, `drop_taxon`, `remove_taxon`,
    `discard_taxon`, `__delitem__`, etc.)

* In most cases, a simple global search-and-replace of "TaxonSet" with
  "TaxonNamespace" and "`taxon_set`" with "`taxon_namespace`" should be
  sufficient to bring existing code into line with DendroPy 4.

* For legacy support, a class called :class:`TaxonSet` exists. This derives with no
  modifications from :class:`TaxonNamespace`. Instantiating objects of this class
  will result in warnings being emitted. As long as usage of :class:`TaxonSet` does
  conforms to the above API change notes, old or legacy code should continue
  to work unchanged (albeit, with some warning noise). This support is
  temporary and will be removed in upcoming releases: code should update to
  using :class:`TaxonNamespace` as soon as expedient.

* For legacy support, "`taxon_set`" continues to be accepted and processed as
  an attribute name and keyword argument synonymous with "`taxon_namespace`".
  Usage of this will result in warnings being emitted, but code should
  continue to function as expected. This support is temporary and will be
  removed in upcoming releases: code should update to using
  "`taxon_namespace`" as soon as expedient.

The :class:`Node` Class
=======================

* Constructor now only accepts keyword arguments (and ``oid`` is *not* one of them!).

The :class:`Edge` Class
=======================

* Constructor now only accepts keyword arguments (and ``oid`` is *not* one of them!).

* Because `tail_node` is no longer an independent attribute but a dynamic
  property, bound to :attr:`Node._parent_node` attribute of the `head_node`
  (see below), the :class:`Edge` constructor does *not* accept ``tail_node`` as
  an argument.

* The `tail_node` of an :class:`Edge` object is now a dynamic property,
  referencing the :attr:`Node._parent_node` attribute of the
  :attr:`Edge._head_node` of the :class:`Edge` object. So, now updating
  :attr:`Edge._tail_node` of an :class:`Edge` object will set the
  :attr:`Node._parent_node` of its :attr:`Edge._head_node` to the new value,
  and vice versa.  This avoids the need for independent book-keeping logic to
  ensure that :attr:`Node._parent_node` and :attr:`Edge._tail_node` are always
  synchornized to reference the same :class:`Node` object and all the potential
  errors this might cause.


The :class:`Tree` Class
=======================

* :meth:`Tree.nodes()` : sorting option removed; use `sorted(tree.nodes())` instead.

* `Tree.node_set()` : removed; use `set(tree.nodes())` instead.

* `Tree.edge_set()` : removed; use `set(tree.edges())` instead.

* For consistency with :meth:`Tree.preorder_node_iter()`,
  :meth:`Tree.postorder_node_iter()`, a number of iteration methods have been renamed.

    +--------------------------------+-------------------------------------+
    | DendroPy 3                     | DendroPy 4                          |
    +--------------------------------+-------------------------------------+
    | `Tree.level_order_node_iter()` | :meth:`Tree.levelorder_node_iter()` |
    +--------------------------------+-------------------------------------+
    | `Tree.level_order_edge_iter()` | :meth:`Tree.levelorder_edge_iter()` |
    +--------------------------------+-------------------------------------+
    | `Node.level_order_iter()`      | :meth:`Node.levelorder_iter()`      |
    +--------------------------------+-------------------------------------+
    | `Edge.level_order_iter()`      | :meth:`Edge.levelorder_iter()`      |
    +--------------------------------+-------------------------------------+
    | `Tree.age_order_node_iter()`   | :meth:`Tree.ageorder_node_iter()`   |
    +--------------------------------+-------------------------------------+
    | `Tree.age_order_edge_iter()`   | :meth:`Tree.ageorder_edge_iter()`   |
    +--------------------------------+-------------------------------------+
    | `Node.age_order_iter()`        | :meth:`Node.ageorder_iter()`        |
    +--------------------------------+-------------------------------------+
    | `Edge.age_order_iter()`        | :meth:`Edge.ageorder_iter()`        |
    +--------------------------------+-------------------------------------+
    | `Tree.leaf_iter()`             | :meth:`Tree.leaf_node_iter()`       |
    +--------------------------------+-------------------------------------+

    The old names are still supported for now (with warnings being emitted),
    but new code should start using the newer names.  In additon, support for
    in-order or infix tree traversal has been added:
    :meth:`Tree.inorder_node_iter`, :meth:`Tree.inorder_edge_iter()`.

