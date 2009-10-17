.. include:: common.rst

********************************
|dendropylogo| DendroPy Cookbook
********************************

.. contents:: Contents


Introduction
============


Reading and Writing Data
========================

Reading Data from a File
------------------------

The ``read()`` method of the ``datasets.Dataset`` object is the primary way of loading data from a file into |DendroPy|_.
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

The ``datasets.Dataset`` object "``d``" will now contain all the data in "primates.tre"---taxa, trees, and characters.
If the data contains only trees or only characters, then the corresponding |DendroPy|_ ``Dataset`` will contain only trees or characters, respectively.
If the data format does not specify an explicit taxon block (e.g., NEWICK, PHYLIP, FASTA, or a technically-invalid but often seen incomplete NEXUS variant), then one will be automatically created and associated with all taxon-linked elements (character blocks, trees blocks, and trees) of the data.

Each ``Dataset`` object has three attributes:

    * ``taxa_blocks`` : a list of ``TaxaBlock`` objects corresponding to the (one or more) taxa blocks in, or implied by, the data
    * ``trees_blocks``: a list of ``TreesBlock`` objects corresponding to the (zero or more) sets of trees in the data
    * ``char_blocks``: a list of ``CharactersBlock`` objects corresponding to the (zero or more) sets of character matrices in th data
    
Most file formats will only result in at most one ``TaxaBlock``, one ``TreesBlock`` and one ``CharactersBlock`` objects in each of the respective lists.
For example, reading a FASTA or PHYLIP file will result in single-element ``taxa_blocks`` and ``char_blocks`` lists.
Similarly, reading a NEWICK file will result in single-element ``taxa_blocks`` and ``trees_blocks`` lists.
A standard NEXUS file will result in a single-element ``taxa_blocks`` list, and either empty or single-element ``char_blocks`` and ``trees_blocks`` lists.

Each ``TaxaBlock`` and ``TreesBlock`` object is, in turn, a specialized list, with ``Taxon`` and ``Tree`` elements respectively.
Each ``Tree`` object consists of nodes (``Node`` objects) and branches (``Edge`` objects), with a ``Taxon`` object from the associated ``TaxaBlock`` assigned to the ``taxon`` attribute of each leaf node.
Internal typically have their ``taxon`` attribute set to ``None``, but this need not neccessarily be the case. Each ``CharactersBlock`` behaves like a dictionary that maps ``Taxon`` objects in its associated ``TaxaBlock`` to vectors of character data.

Almost every object has a ``label`` attribute, which is a plain |Python|_ string. 
It is important to distinguish between the string label of an object and the object itself. 
For example, a NEXUS file may contain a tree which includes a taxon label "Agkistrodon".
When this file is read by |DendroPy|_, a ``Taxon`` object will be created with its ``label`` attribute set to "Agkistrodon".
However, the ``Tree`` object created will not have the ``label`` attribute of the corresponding node set to anything.
Instead, the ``taxon`` attribute of the node will point to a ``Taxon`` object with the label "Agkistrodon" in its associated ``TaxaBlock``.
Some file formats, such as NEXUS or NeXML, allow for labels to be associated with nodes independentally of taxa.
For example, the NEXUS specification allows for internal node labels.
*These* labels *will* result in the ``label`` attribute being set on the corresponding nodes of the DendroPy ``Tree`` object.


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
         
The following script performs something I find *very* useful: it reads in a FASTA file downloaded from GenBank, and writes out the data in NEXUS format, transforming the highly-informative but also verbose GenBank labels to something that is meaningful and yet valid for direct use in most phylogenetic programs::

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

Once a ``Dataset`` object has been instantiated, by examining the lengths of the lists of ``taxa.TaxaBlock``, ``trees.TreeBlock`` and ``characters.CharactersBlock`` objects we can determined how many of each kind are there::

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




