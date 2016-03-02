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
    - |StandardCharacterMatrix|, for discrete-value data
    - |ContinuousCharacterMatrix|, for continuous-valued data

The |ContinuousCharacterMatrix| class represents its character values directly.
Typically, all other classes represent its character values as special :class:`~dendropy.datamodel.charstatemodel.StateIdentity` instances, *not* as strings.
So, for example, the DNA character "A" is modeled by a special :class:`~dendropy.datamodel.charstatemodel.StateIdentity` instance (created by the DendroPy library).
While it is represented by the string "A", and can be converted to the string and back again, it is not the same as the string "A".
Each discrete |CharacterMatrix| instance has one or more :class:`~dendropy.datamodel.charstatemodel.StateAlphabet` instances associated with it that manage the collection of letters that make up the character data.
In the case of, e.g. DNA, RNA, protein and other specialized discrete data, this are pre-defined by DendroPy: ``dendropy.DNA_STATE_ALPHABET``, ``dendropy.RNA_STATE_ALPHABET``, etc.
In the case of "standard" character data, these are created for each matrix separately. Facilities are provided for the creation of custom state alphabets and for the sharing of state alphabets between different |StandardCharacterMatrix| instances.

Reading and Writing Character Data
==================================

As with most other phylogenetic data objects, objects of the |CharacterMatrix|-derived classes support the "|get|" factory method to populate objects from a data source.
This method takes a data source as the first keyword argument and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` ("``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``", etc.) as the second::

    import dendropy
    dna1 = dendropy.DnaCharacterMatrix.get(file=open("pythonidae.fasta"), schema="fasta")
    dna2 = dendropy.DnaCharacterMatrix.get(url="http://purl.org/phylo/treebase/phylows/matrix/TB2:M2610?format=nexus", schema="nexus")
    aa1 = dendropy.ProteinCharacterMatrix.get(file=open("pythonidae.dat"), schema="phylip")
    std1 = dendropy.StandardCharacterMatrix.get(path="python_morph.nex", schema="nexus")
    std2 = dendropy.StandardCharacterMatrix.get(data=">t1\n01011\n\n>t2\n11100", schema="fasta")


The "|write|" method allows you to write the data of a |CharacterMatrix| to a file-like object or a file path::

    dna1 = dendropy.DnaCharacterMatrix.get(file=open("pythonidae.nex"), schema="nexus")
    dna1.write(path="out.nexml", schema="nexml")
    dna1.write(file=open("out.fasta", schema="fasta")

You can also represent the data as a string using the :meth:`as_string` method::

    dna1 = dendropy.DnaCharacterMatrix.get(file=open("pythonidae.nex"), schema="nexus")
    s = dna1.as_string(schema="fasta")
    print(s)


In addition, fine-grained control over the reading and writing of data is available through various keyword arguments as described in the :doc:`/primer/reading_and_writing` section.

Creating a Character Data Matrix from a Dictionary of Strings
=============================================================

The :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.from_dict` factory method creates a new |CharacterMatrix| from a dictionary mapping taxon labels to sequences represented as strings::

    import dendropy
    d = {
            "s1" : "TCCAA",
            "s2" : "TGCAA",
            "s3" : "TG-AA",
    }
    dna = dendropy.DnaCharacterMatrix.from_dict(d)

Taxon Management with Character Matrices
========================================

Taxon management with |CharacterMatrix|-derived objects work very much the same as it does with |Tree| or |TreeList| objects every time a |CharacterMatrix|-derived object is independentally created or read, a new |TaxonNamespace| is created, unless an existing one is specified.
Thus, again, if you are creating multiple character matrices that refer to the same set of taxa, you will want to make sure to pass each of them a common |TaxonNamespace| reference::

    import dendropy
    taxa = dendropy.TaxonNamespace()
    dna1 = dendropy.DnaCharacterMatrix.get(
        path="pythonidae_cytb.fasta",
        schema="fasta",
        taxon_namespace=taxa)
    prot1 = dendropy.ProteinCharacterMatrix.get(
        path="pythonidae_morph.nex",
        schema="nexus",
        taxon_namespace=taxa)
    trees = dendropy.TreeList.get(
        path="pythonidae.trees.nex",
        schema="nexus",
        taxon_namespace=taxa)

Concatenating Multiple Data Matrices
====================================

A new |CharacterMatrix| can be created from multiple existing matrices using the :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatentate()` factory method, which takes a list or an iterable of |CharacterMatrix| instances as an argument.
        All the CharacterMatrix objects in the list must be of the
        same type, and share the same TaxonNamespace reference. All taxa
        must be present in all alignments, all all alignments must
        be of the same length. Component parts will be recorded as
        character subsets.

For example:

.. literalinclude:: /examples/char_mat_concat.py


results in ::

    d1: 12 sequences, 231 characters
    d2: 12 sequences, 231 characters
    d3: 12 sequences, 231 characters
    d_all: 12 sequences, 693 characters
    Subsets: {'locus002': <dendropy.datamodel.charmatrixmodel.CharacterSubset object at 0x101d792d0>, 'locus000': <dendropy.datamodel.charmatrixmodel.CharacterSubset object at 0x101d79250>, 'locus001': <dendropy.datamodel.charmatrixmodel.CharacterSubset object at 0x101d79290>}

You can instantiate a concatenated matrix from multiple sources using the :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatentate_from_paths()` or :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatentate_from_streams()` factory methods:

.. literalinclude:: /examples/char_mat_concat2.py

Sequence Management
===================

A range of methods also exist for importing data from another matrix object.
These vary depending on how "new" and "existing" are treated.  A "new"
sequence is a sequence in the other matrix associated with a |Taxon|
object for which there is no sequence defined in the current matrix.  An
"existing" sequence is a sequence in the other matrix associated with a
|Taxon| object for which there *is* a sequence defined in the
current matrix.

+---------------------------------+---------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
|                                 | New Sequences: IGNORED                                                          | New Sequences: ADDED                                                           |
+=================================+=================================================================================+================================================================================+
| Existing Sequences: IGNORED     | [NO-OP]                                                                         | :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.add_sequences()`    |
+---------------------------------+---------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| Existing Sequences: OVERWRITTEN | :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.replace_sequences()` | :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.update_sequences()` |
+---------------------------------+---------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| Existing Sequences: EXTENDED    | :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.extend_sequences()`  | :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.extend_matrix()`    |
+---------------------------------+---------------------------------------------------------------------------------+--------------------------------------------------------------------------------+

More information cane be found in the source documentation:

-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.add_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.replace_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.update_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.extend_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.extend_matrix()`

In addition there are methods for selecting removing sequences:

-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.remove_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.discard_sequences()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.keep_sequences()`

As well as "filling out" a matrix by adding columns or rows:

-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.fill_taxa()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.fill()`
-   :meth:`~dendropy.datamodel.charmatrixmodel.CharacterMatrix.pack()`

Accessing Data
==============

A |CharacterMatrix| behaves very much like a dictionary, where the "keys" are |Taxon| instances, which can be dereferenced using the instance itself, the taxon label, or the index of the taxon in the collection (note: this is *not* neccessarily the same as the accession index, which is the basis for bipartition collection).

For example:

.. literalinclude:: /examples/chars_access1.py

You can also iterate over the matrix in a number of ways:

.. literalinclude:: /examples/chars_access2.py


The "values" return by dereferencing the "keys" of a |CharacterMatrix|  objects are |CharacterDataSequence| objects.
Objects of this class behave very much like lists, where the elements are either numeric values for |ContinuousCharacterMatrix| matrices:

.. literalinclude:: /examples/chars_access3.py

or |StateIdentity| instances for all other types of matrices:

.. literalinclude:: /examples/chars_access4.py

As can be seen, you can use :meth:`~dendropy.datamodel.charmatrixmodel.CharacterDataSequence.values()` to get a list of the values of the sequence directly, :meth:`~dendropy.datamodel.charmatrixmodel.CharacterDataSequence.symbols_as_list()` to get a list of the values represented as strings, and :meth:`~dendropy.datamodel.charmatrixmodel.CharacterDataSequence.symbols_as_string()` to get the string representation of the whole sequence.
