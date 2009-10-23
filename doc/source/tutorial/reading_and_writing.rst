************************
Reading and Writing Data
************************

Reading Data from a File
========================

The ``read()`` method of the ``datasets.Dataset`` object is the primary way of loading data from a file into |DendroPy|_.
It takes two parameters, a file handle and a (case-insensitive) string specifying a format, which can be one of:

    * "NEXUS"
    * "NEWICK"
    * "NEXML"
    * "PHYLIP"
    * "DNAFASTA"
    * "RNAFASTA"

For example, to load the data in the NEXUS-format file, "primates.tre":

.. topic:: Reading in a Data File
    :class: code-recipe
    
    :: 

        >>> from dendropy import datasets
        >>> d = datasets.Dataset()
        >>> d.read( open("primates.tre", "rU"), "NEXUS" )
        <dendropy.datasets.Dataset object at 0x2a26f0>

The ``datasets.Dataset`` object "``d``" will now contain all the data in "primates.tre"---taxa, trees, and characters.

Reading Data from Strings
=========================

You can also create datasets directly from string, using the :meth:`from_string` method::

    >>> from dendropy import datasets
    >>> d = dataset.Dataset()
    >>> d.from_string("(a,(b,c)); ((a,b),c);", format="newick")
        <dendropy.datasets.Dataset object at 0x755250>

:obj:`d.taxa_blocks` now will contain a single :class:`TaxaBlock` element, with three :class:`Taxon` object with labels "a", "b" and "c", while :obj:`d.trees_blocks` will have a single :class:`TreesBlock` element, which will have two :class:`Tree` objects.

Writing Data to a File
=======================
The ``write()`` method of the ``Dataset`` object writes the data to a file. As with ``read()``, it takes two arguments: a file handle and a case-insenstive string specifying the format.
The following writes out the data we just read into a NEWICK format file:

.. topic:: Writing to a Data File
    :class: code-recipe
    
    :: 

        >>> d.write(open("primates.newick.tre", "w"), "NEWICK")

Translating Between Different Formats
=====================================
Actual serialization formats are (largely) opaque to the DendroPy data model, with all deserialization/serialization handled by specialized dataset readers and writers respectively.
So reading and writing in different formats is simply a matter of calling ``read()`` and ``write()`` with the appropriate format specifications, as the examples below demonstrate.

.. topic:: Converting Data from FASTA to NEXUS Format
    :class: code-recipe
    
    :: 

        >>> from dendropy import datasets
        >>> d = datasets.Dataset()
        >>> d.read(open("rana.fasta", "rU"), "DNAFASTA")
        >>> d.write(open("rana.nex", "w"), "NEXUS")
    
|    
    
.. topic:: Converting Data from NEXUS to PHYLIP Format
    :class: code-recipe
    
    :: 
    
        >>> from dendropy import datasets
        >>> d = datasets.Dataset()
        >>> d.read(open("rana.nex", "rU"), "NEXUS")
        >>> d.write(open("rana.dat", "w"), "PHYLIP")
    
|

.. topic:: Converting Data from PHYLIP to FASTA Format
    :class: code-recipe
    
    ::
    
        >>> from dendropy import datasets
        >>> d = datasets.Dataset()
        >>> d.read(open("rana.dat", "rU"), "PHYLIP")
        >>> d.write(open("rana2.fasta", "w"), "FASTA")
             
                  
The following script performs something I find *very* useful: it reads in a FASTA file downloaded from GenBank, and writes out the data in NEXUS format, transforming the highly-informative but also verbose GenBank labels to something that is meaningful and yet valid for direct use in most phylogenetic programs.

.. topic:: Fixing Labels in a GenBank FASTA File
    :class: code-recipe
    
    ::        

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
