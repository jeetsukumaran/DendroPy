*******************************
Working with Character Matrices
*******************************

The |CharacterMatrix| object represents character data in DendroPy.
In most cases, you will not deal with objects of the |CharacterMatrix| class directly, but rather with objects of one of the classes specialized to handle specific data types:

    - |DnaCharacterMatrix|, for DNA nucleotide sequence data
    - |RnaCharacterMatrix|, for RNA nucleodtide sequence data
    - |ProteinCharacterMatrix|, for amino acid sequence data
    - |ContinuousCharacterMatrix|, for continuous-valued data
    - |StandardCharacterMatrix|, for discrete-value data

As with most other phylogenetic data objects, objects of the |CharacterMatrix|-derived classes support the :meth:`get_from_*()` factory and :meth:`read_from_*()` instance methods to populate objects from a data source.
The following examples simultaneously instantiate and populate |CharacterMatrix| objects of the appropriate type from various file data sources::

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_cytb.nex', 'nexus')
    >>> rna = dendropy.DnaCharacterMatrix.get_from_path('hiv1_env.nex', 'nexus')
    >>> aa = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_mos.nex', 'nexus')
    >>> cv = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_sizes.nex', 'nexus')
    >>> sm = dendropy.DnaCharacterMatrix.get_from_path('pythonidae_skull.nex', 'nexus')

The :meth:`read_from_*()` instance methods **replace** the calling object with data from the data source, overwriting existing data::

    >>> import dendropy
    >>> dna = dendropy.DnaCharacterMatrix()
    >>> dna.read_from_path('pythonidae_cytb.nex', 'nexus')
    >>> dna.read_from_path('pythonidae_rag1.nex', 'nexus')

The second :meth:`read_from_*()` will result in the ``dna`` object being re-populated with data from the file ``pythonidae_rag1.nex``.

.. _Customizing_Character_Creation_and_Reading:

Customizing Character Matrix Creation and Reading
=================================================

Using a Specific |TaxonSet|
---------------------------
Passing a |TaxonSet| object using the ``taxon_set`` argument when using the meth:`get_from_*()` or :meth:`read_from_*()` methods results in the |CharacterMatrix| object being bound to the specified |TaxonSet| object.

Custom Handling of Underscores
------------------------------
With NEXUS and NEWICK data sources, you can also specify ``preserve_underscores=True``.
The NEXUS standard dictates that underscores are equivalent to spaces, and thus any underscore found in any unquoted label in a NEXUS/NEWICK data source will be substituted for spaces.
Specifying ``preserve_underscores=True`` will force DendroPy to keep the underscores. More details on using this keyword to manage taxon references and mapping can be found in here: :ref:`Taxon_Label_Mapping`.

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
