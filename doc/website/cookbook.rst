.. |dendropylogo| image:: dendropy-logo.png
                  :class: dendropylogo
    
********************************
|dendropylogo| DendroPy Cookbook
********************************

.. contents:: Contents


Introduction
============

The DendroPy Data Model
=======================

Dataset Objects
---------------

The DendroPy class ``datasets.Dataset`` is the primary dataset class.
All phylogenetic data, including taxa, trees, and characters, are managed by objects of this class.
A ``Dataset`` object has three attributes:

    * "``taxa_blocks``": a list of ``taxa.TaxaBlock`` objects
    * "``trees_blocks``": a list of ``trees.TreeBlock`` objects
    * "``char_blocks``": a list of ``characters.CharBlock`` objects
    
Each ``taxa.TaxaBlock`` object is a specialized list of ``taxa.Taxon`` objects, while each ``trees.TreeBlock`` is similarly a specialized list of ``trees.Tree`` objects.
Each ``characters.CharBlock``, on the other hand, can be seen as a specialized dictionary which maps each ``taxa.Taxon`` to a vector of characters representing a character data sequence associated with that taxon.

The DendroPy data model can be considered a taxon-driven data model, in that that all tree and character data are partitioned by taxa.
That is, each and every ``TreeBlock`` and ``CharBlock`` are associated with one, and only one, ``TaxaBlock`` object, but each ``TaxaBlock`` object may be associated with more than one (or none) ``TreeBlock`` and ``CharBlock`` objects.

This taxon-driven partition applies to the individual ``TreeBlock`` objects as well as ``CharBlock`` objects: every ``TreeBlock`` and ``CharBlock`` object is associated with a ``TaxaBlock`` object that determines the full set of taxa represented in the tree and character data respectively.

Adding a ``Tree`` object to a ``TreeBlock`` list automatically normalizes the ``TaxaBlock`` associated with the ``Tree`` object: the ``TaxaBlock`` associated with the ``TreeBlock`` list will be assigned to the ``Tree`` object, and ``taxa.Taxon`` objects on the ``Tree`` object will be replaced with corresponding ``taxa.Taxon`` objects of the ``TreeBlock`` list, with identity established by the taxon label. If there are ``taxa.Taxon`` objects in the ``Tree`` object with labels that do not match with the labels of any ``Taxon`` objects in the ``TaxaBlock`` of the ``TreeBlock`` list, then a new ``Taxon`` object with the label will be created and added to the ``TaxaBlock`` of the ``TreeBlock`` list.

Tree Objects
-------------

DendroPy ``Tree`` objects are rich objects. 
Every ``Tree`` object has a ``seed_node`` attribute, which points to the starting node of the tree.
This is the rooted if the tree is rooted (``Tree.is_rooted=True``), but is an artifical node if it is not.
Every node on the ``Tree`` object (including ``seed_node``) is an object of the ``Trees.Node`` class.
Each ``Trees.Node`` has, among others, the following attributes:

    * ``edge`` : a ``Trees.Edge`` object representing the edge subtending this node
    * ``parent_node`` : a ``Node`` object that is the parent of the current node (``None`` for the root or ``seed_node``)
    * ``label`` : a text label associated with this node
    * ``taxon`` : a ``Taxon`` object representing the taxon associated with this node (typically, internal nodes will have this set to ``None``, but this is not always the case)

Reading and Writing Data
========================

Reading Data from a File
------------------------

The ``read()`` method of the ``datasets.Dataset`` object is the primary way of loading data from a file to a dataset.
It takes two parameters, a file handle and a (case-insensitive) string specifying a format, which can be one of:

    * "NEXUS"
    * "NEWICK"
    * "NEXML"
    * "PHYLIP"
    * "DNAFASTA"
    * "RNAFASTA"    

For example, to load the data in the NEXUS-format file, "primates.tre"::

    >>> from dendropy import datasets
    >>> d = datasets.Dataset()
    >>> d.read( open("primates.tre", "rU"), "NEXUS" )
    <dendropy.datasets.Dataset object at 0x2a26f0>

The ``datasets.Dataset`` object "``d``" will now contain all the data in "primates.tre".

Writing Data to a File
-----------------------
The ``write()`` method of the ``Dataset`` object writes the data to a file. As with ``read()``, it takes two arguments: a file handle and a case-insenstive string specifying the format.
The following writes out the data we just read into a NEWICK format file::

    >>> d.write(open("primates.newick.tre", "w"), "NEWICK")

Translating Between Different Formats
-------------------------------------
Actual serialization formats are (largely) opaque to the DendroPy data model, with all deserialization/serialization handled by specialized dataset readers and writers respectively.
So reading and writing in different formats is simply a matter of calling ``read()`` and ``write()`` with the appropriate format specifications, as the examples below demonstrate.

FASTA to NEXUS::

    >>> from dendropy import datasets
    >>> d = datasets.Dataset()
    >>> d.read(open("rana.fasta", "rU"), "DNAFASTA")
    >>> d.write(open("rana.nex", "w"), "NEXUS")
    
NEXUS to PHYLIP::

    >>> from dendropy import datasets
    >>> d = datasets.Dataset()
    >>> d.read(open("rana.nex", "rU"), "NEXUS")
    >>> d.write(open("rana.dat", "w"), "PHYLIP")
    
PHYLIP to FASTA::

    >>> from dendropy import datasets
    >>> d = datasets.Dataset()
    >>> d.read(open("rana.dat", "rU"), "PHYLIP")
    >>> d.write(open("rana2.fasta", "w"), "DNAFASTA")
         
The following script performs something I find *very* useful: it reads in a FASTA file and writes out the data in NEXUS format, transforming the labels to something that is meaningful and yet valid::

    #! /usr/bin/env python
    
    import re
    import sys
    from dendropy import datasets
    
    fd = datasets.Dataset()
    fd.read(open("python_cytb.fasta", "rU"), "DNAFASTA")
    pattern = re.compile("gi\|.+\|.+\|(.+)\|\S* ([\w\.]+) ([\w\.]+) (\w+).*")
    for t in fd.taxa_blocks[0]:
        m = pattern.match(t.label)
        t.label = m.groups(1)[0] + "_" + m.groups(1)[1] + "_" + m.groups(1)[2]
        sys.stderr.write(t.label + "\n")
    fd.write(open("python_cytb.nexus", "w"), "NEXUS")


Accessing and Querying Data
============================

Once a ``Dataset`` object has been instantiated, by examining the lengths of the lists of ``taxa.TaxaBlock``, ``trees.TreeBlock`` and ``characters.CharBlock`` objects we can determined how many of each kind are there::

    >>> from dendropy import datasets
    >>> d = datasets.Dataset()
    >>> d.read( open("primates.tre", "rU"), "NEXUS" )
    <dendropy.datasets.Dataset object at 0x2a26f0>
    >>> len(d.taxa_blocks)
    1
    >>> len(d.trees_blocks)
    1
    >>> len(d.char_blocks)
    0    

The first, and only, element in the list of taxa blocks is a ``TaxaBlock`` object, which is in turn a specialized list that contains all the taxa in the file::

    >>> d.taxa_blocks[0]
    [<DendroPy Taxon: 'Lemur catta'>, <DendroPy Taxon: 'Homo sapiens'>, <DendroPy Taxon: 'Pan'>, <DendroPy Taxon: 'Gorilla'>, <DendroPy Taxon: 'Pongo'>, <DendroPy Taxon: 'Hylobates'>, <DendroPy Taxon: 'Macaca fuscata'>, <DendroPy Taxon: 'Macaca mulatta'>, <DendroPy Taxon: 'Macaca fascicularis'>, <DendroPy Taxon: 'Macaca sylvanus'>, <DendroPy Taxon: 'Saimiri sciureus'>, <DendroPy Taxon: 'Tarsius syrichta'>]
    
And similarly for the ``trees_blocks`` attribute of the dataset::

    >>> d.trees_blocks[0]
    [<dendropy.trees.Tree object at 0x5a9690>, <dendropy.trees.Tree object at 0x5a9730>]

Iterating Through Taxa
----------------------
The following snippet loops over the taxa in the first taxa block, printing their labels::

    >>> for t in d.taxa_blocks[0]:
    ...     print(t.label)
    ... 
    Lemur catta
    Homo sapiens
    Pan
    Gorilla
    Pongo
    Hylobates
    Macaca fuscata
    Macaca mulatta
    Macaca fascicularis
    Macaca sylvanus
    Saimiri sciureus
    Tarsius syrichta

Iterating Through Trees
-----------------------
The same approach works for the trees::

    >>> for t in d.trees_blocks[0]:
    ...     print(t.label)
    ... 
    rep.1
    rep.1000

We can also inspect the NEWICK string representations of the trees::

    >>> for t in d.trees_blocks[0]:
    ...     print(t.compose_newick())
    ... 
    ((((('Macaca fascicularis':0.1,'Tarsius syrichta':0.1):0.1,'Saimiri sciureus':0.121635):0.089589,(('Macaca fuscata':0.1,Gorilla:0.1):0.1,(('Macaca sylvanus':0.1,Pan:0.1):0.1,Hylobates:0.1):0.1):0.100676):0.1,'Homo sapiens':0.1):0.1,('Macaca mulatta':0.1,Pongo:0.1):0.1,'Lemur catta':0.1)
    ('Tarsius syrichta':0.247169,(('Saimiri sciureus':0.325537,(('Macaca fascicularis':0.065018,('Macaca mulatta':0.022964,'Macaca fuscata':0.020959):0.02792):0.028642,'Macaca sylvanus':0.088559):0.246816):0.019503,((Pongo:0.093129,(('Homo sapiens':0.044705,Pan:0.082301):0.011332,Gorilla:0.061149):0.066643):0.068598,Hylobates:0.154276):0.090646):0.243449,'Lemur catta':0.258383)

Tree Traversal
--------------

Trees can be traversed in pre-order, post-order, or level-order, over nodes or edges.




