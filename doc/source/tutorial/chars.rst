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
These methods take a file-like object (in the case of :meth:`write_to_stream()`) or a string specifying a filepath (in the case of :meth:`write_to_path()`) as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the second argument.

The following example reads a FASTA-formatted file and writes it out to a a NEXUS-formatted file:

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> dna.write_to_path('pythonidae_cytb.nexus', 'nexus')

Fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

Composing a String
------------------

If you do not want to actually write to a file, but instead simply need a string representing the data in a particular format, you can call the instance method :meth:`as_string()`, passing a :ref:`schema specification string <Specifying_the_Data_Writing_Format>` as the first argument:

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.fasta', 'dnafasta')
    >>> s = dna.as_string('nexus')
    >>> print(s)

As above, fine-grained control over the output format can be specified using :ref:`keyword arguments <Customizing_the_Data_Writing_Format>`.

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

