******************
Character Matrices
******************

Types of Character Matrices
===========================

The |CharacterMatrix| object represents character data in DendroPy.
In most cases, you will not deal with objects of the |CharacterMatrix| class directly, but rather with objects of one of the classes specialized to handle specific data types:

    - |DnaCharacterMatrix|, for DNA nucleotide sequence data
    - |RnaCharacterMatrix|, for RNA nucleodtide sequence data
    - |ProteinCharacterMatrix|, for amino acid sequence data
    - |ContinuousCharacterMatrix|, for continuous-valued data
    - |StandardCharacterMatrix|, for discrete-value data

|CharacterMatrix| Creating and Reading
======================================

As with most other phylogenetic data objects, objects of the |CharacterMatrix|-derived classes support the :meth:`get_from_*()` factory and :meth:`read_from_*()` instance methods to populate objects from a data source.
These methods take a data source as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``", etc.) as the second, as well as optional :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>` to customize the reading behavior.

Creating a New |CharacterMatrix| from a Data Source
---------------------------------------------------

The following examples simultaneously instantiate and populate |CharacterMatrix| objects of the appropriate type from various file data sources::

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.nex', 'nexus')
    >>> rna = dendropy.DnaCharacterMatrix.get_from_path('hiv1_env.nex', 'nexus')
    >>> aa = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_mos.nex', 'nexus')
    >>> cv = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_sizes.nex', 'nexus')
    >>> sm = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_skull.nex', 'nexus')

Repopulating a |CharacterMatrix| from a DataSource
--------------------------------------------------

The :meth:`read_from_*()` instance methods **replace** the calling object with data from the data source, overwriting existing data::

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix()
    >>> dna.read_from_path('pythonidae_cytb.nex', 'nexus')
    >>> dna.read_from_path('pythonidae_rag1.nex', 'nexus')

The second :meth:`read_from_*()` will result in the ``dna`` object being re-populated with data from the file ``pythonidae_rag1.nex``.

|CharacterMatrix| Saving and Writing
====================================

Writing to Files
----------------

The :meth:`write_to_stream()`, and :meth:`write_to_path()` instance methods allow you to write the data of a |CharacterMatrix| to a file-like object or a file path respectively.
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a format or schema specification string as the second argument.

The following example reads a FASTA-formatted file and writes it out to a a NEXUS-formatted file:

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> dna.write_to_path('pythonidae_cytb.nexus', 'nexus')

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a schema or format specification string as the first argument:

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> s = dna.as_string('nexus')
    >>> print(s)

Taxon Management with Character Matrices
========================================

Taxon management with |CharacterMatrix|-derived objects work very much the same as it does with |Tree| or :class:`~dendropy.dataobject.tree.TreeList objects`: every time a |CharacterMatrix|-derived object is independentally created or read, a new |TaxonSet| is created, unless an existing one is specified.
Thus, again, if you are creating multiple character matrices that refer to the same set of taxa, you will want to make sure to pass each of them a common |TaxonSet| reference::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> dna1 = dendropy.DnaCharacterMatrix.get_from_path("pythonidae_cytb.fasta", "dnafasta", taxon_set=taxa)
    >>> std1 = dendropy.ProteinCharacterMatrix.get_from_path("pythonidae_morph.nex", "nexus", taxon_set=taxa)


Accessing Data
==============
Each sequence for a particular |Taxon| object is organized into a |CharacterDataVector| object, which, in turn, is a list of |CharacterDataCell| objects.
You can retrieve the |CharacterDataVector| for a particular taxon by passing the corresponding |Taxon| object, its label, or its index to the |CharacterMatrix| object.
Thus, to get the character sequence vector associated with the first taxon ("``Python regius``") from the data source ``pythonidae_cytb.fasta``:

    >>> from dendropy import DnaCharacterMatrix
    >>> cytb = DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> v1 = cytb[0]
    >>> v2 = cytb['Python regius']
    >>> v3 = cytb[cytb.taxon_set[0]]
    >>> v1 == v2 == v3
    True

Population Genetic Summary Statistics
=====================================

The :mod:`popgenstat` module provides functions that calculate some common population genetic summary statistics.

For example, given a |DnaCharacterMatrix| as an argument, the :func:`~dendropy.popgenstat.num_segregating_sites()` function returns the raw number of segregating sites, :func:`~dendropy.popgenstat.average_number_of_pairwise_differences()` returns the average number of pairwise differences, and :func:`~dendropy.popgenstat.nucleotide_diversity()` returns the nucleotide diversity.

More complex statistics are provided by the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` class.
Objects of this class are instantatiated with two lists of |CharacterDataVector| objects as arguments, each representing a sample of DNA sequences drawn from two distinct but related populations.
Once instantiated, the following attributes of the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` object are available:

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences`
            The average number of pairwise differences between every sequence across both populations.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_between`
            The average number of pairwise differences between every sequence between both populations.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_within`
            The average number of pairwise differences between every sequence within each population.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_net`
            The net number of pairwise differences.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.num_segregating_sites`
            The number of segregating sites.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.wattersons_theta`
            Watterson's theta.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.wakeleys_psi`
            Wakeley's psi.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.tajimas_d`
            Tajima's D.

The following example calculates the suite of population genetic summary statistics for sequences drawn from two populations of sticklebacks.
The original data consists of 23 sequences, with individuals from Eastern Pacific populations identified by their taxon labels beginning with "``EPAC``" and individuals from Western Pacific populations identified by their taxon labels beginning with "``WPAC``".
The taxon labels thus are used as the basis for sorting the sequences into the required lists of |CharacterDataVector| objects, ``p1`` and ``p2``.

.. literalinclude:: /examples/pgstats1.py
    :linenos:

Lines 6-12 build up the two lists of |CharacterDataVector| objects by sorting the original sequences into their source populations based on the taxon label (with operational taxonomic units with labels beginning with "``EPAC``" coming from the Eastern Pacific, and assigned to the list ``p1``, while those that begin with "``WPAC``" coming from the Western Pacific, and assigned to the list ``p2``).
These lists are then passed as the instantiation arguments to the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` constructor in line 14.
The calculations are performed immediately, and the results reported in the following lines.
