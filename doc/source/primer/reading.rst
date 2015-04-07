*************************
Reading Phylogenetic Data
*************************

Creating and Populating New Objects
===================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes all
support a "|get|" factory class-method that instantiates an object of the given
class from a data source.
This method takes, at a minumum, two keyword arguments that specify the
*source* of the data and the *schema* (or format) of the data.

The source can be:

    -   a path to a file (specified using the keyword argument "``path``")
    -   a file or a file-like object opened for reading (specified using the keyword argument ``"file"``)
    -   a string value giving the data directly (specified using the keyword argument ``"value"``)
    -   or a URL (specified using the keyword argument ``"url"``)

The schema is specified using the keyword argument ``"schema"``, and takes a string value that identifies the format of data.
This ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" can be one of: "``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``".
Not all formats are supported for reading, and not all formats make sense for particular objects (for example, it would not make sense to try and instantiate a |Tree| or |TreeList| object from a FASTA-formatted data source).

.. A ":term:`schema`" is DendroPy-speak for "format" (we cannot use the argument or variable name "format" for this in library, because this is a Python built-in, and hence we use "schema" and adopted this terminology for consistency), and is specified using one of a set of predefined string values.

For example:

.. code-block:: python

    import dendropy

    tree1 = dendropy.Tree.get(path="mle.tre", schema="newick")
    tree2 = dendropy.Tree.get(file=open("mle.nex", "r"), schema="nexus")
    tree3 = dendropy.Tree.get(value="((A,B),(C,D));", schema="newick")
    tree4 = dendropy.Tree.get(url="http://api.opentreeoflife.org/v2/study/pg_1144/tree/tree2324.nex", schema="nexus")

    tree_list1 = dendropy.TreeList.get(path="pythonidae.mcmc.nex", schema="nexus")
    tree_list2 = dendropy.TreeList.get(file=open("pythonidae.mcmc.nex", "r"), schema="nexus")
    tree_list3 = dendropy.TreeList.get(value="(A,(B,C));((A,B),C);", "r"), schema="newick")

    dna1 = dendropy.DnaCharacterMatrix.get(file=open("pythonidae.fasta"), schema="fasta")
    dna2 = dendropy.DnaCharacterMatrix.get(url="http://purl.org/phylo/treebase/phylows/matrix/TB2:M2610?format=nexus", schema="nexus")
    aa1 = dendropy.ProteinCharacterMatrix.get(file=open("pythonidae.dat"), schema="phylip")
    std1 = dendropy.StandardCharacterMatrix.get(path="python_morph.nex", schema="nexus")
    std2 = dendropy.StandardCharacterMatrix.get(value=">t1\n01011\n\n>t2\n11100", schema="fasta")

    dataset1 = dendropy.DataSet.get(path="pythonidae.chars_and_trees.nex", schema="nexus")
    dataset2 = dendropy.DataSet.get(url="http://purl.org/phylo/treebase/phylows/study/TB2:S1925?format=nexml", schema="nexml")

The "|get|" method takes a number of other optional keyword arguments that provide control over how the data is interpreted and processed.
Some of these are specific to the given class (e.g. the "``edge_length_type``" argument when reading tree data into a |Tree|, |TreeList|, or |DataSet| object), while others are specific to the schema or format (e.g., the "``preserve_underscores``" argument when reading |Newick| or |Nexus| data).


.. note::

    These classes also support a "|get_from_methods|" family of factory
    class-methods that can be seen as specializations of the "|get|" method for
    various types of sources (in fact, the "|get|" method is actually a
    dispatcher that calls on one of these methods below for implementation of the functionality):

        :meth:`get_from_stream(src, schema, \*\*kwargs)`
            Takes a file or file-like object opened for reading the data source as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.

        :meth:`get_from_path(src, schema, \*\*kwargs)`
            Takes a string specifying the path to the the data source file as the first argument, and a string specifying the :term:`schema` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.

        :meth:`get_from_string(src, schema, \*\*kwargs)`
            Takes a string containing the source data as the first argument, and a string specifying the :term:`schema` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.

        :meth:`get_from_url(src, schema, \*\*kwargs)`
            Takes a string containing the URL of the data as the first argument, and a string specifying the :term:`schema` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.

    As with the "|get|" method, the additional keywords are specific to the given class or schema type.

Reading and Populating (or Repopulating) Existing Objects
=========================================================

In addition to the "|get_from_methods|" class factory methods, the collection classes (|TreeList|, |TreeArray| and |DataSet|) all support a suite of "|read_from_methods|" *instance* methods that *add* data from external sources to an existing object (as opposed to creating and returning a new object based on an external data source).
These "|read_from_methods|" instance methods have signatures that parallel the "|get_from_methods|" factory methods described above:

    :meth:`read_from_stream(src, schema, **kwargs)`
        Takes a file or file-like object opened for reading the data source as the first argument, and a string specifying the :term:`schema` as the second.

    :meth:`read_from_path(src, schema, **kwargs)`
        Takes a string specifying the path to the the data source file as the first argument, and a string specifying the :term:`schema` as the second.

    :meth:`read_from_string(src, schema, **kwargs)`
        Takes a string specifying containing the source data as the first argument, and a string specifying the :term:`schema` as the second.

When called on an existing |TreeList|, |TreeArray| or |DataSet| objects, these methods *add* the data from the data source to the object.
As with the "|get_from_methods|" methods, the :ref:`schema specification string <Specifying_the_Data_Source_Format>` can be any supported and type-apppropriate :term:`schema`, such as "``nexus``", "``newick``", "``nexml``", "``fasta``", "``phylip``", etc.

For example, the following accumulates post-burn-in trees from several different files into a single |TreeList| object::

    >>> import dendropy
    >>> post_trees = dendropy.TreeList()
    >>> post_trees.read_from_path("pythonidae.nex.run1.t", "nexus", tree_offset=200)
    >>> print(post_trees.description())
    TreeList object at 0x550990 (TreeList5573008): 801 Trees
    >>> post_trees.read_from_path("pythonidae.nex.run2.t", "nexus", tree_offset=200)
    >>> print(post_trees.description())
    TreeList object at 0x550990 (TreeList5573008): 1602 Trees
    >>> post_trees.read_from_path("pythonidae.nex.run3.t", "nexus", tree_offset=200)
    >>> print(post_trees.description())
    TreeList object at 0x550990 (TreeList5573008): 2403 Trees
    >>> post_trees.read_from_path("pythonidae.nex.run4.t", "nexus", tree_offset=200)
    >>> print(post_trees.description())
    TreeList object at 0x5508a0 (TreeList5572768): 3204 Trees

The |TreeList| object automatically handles taxon management, and ensures that all appended |Tree| objects share the same |TaxonNamespace| reference. Thus all the |Tree| objects created and aggregated from the data sources in the example will all share the same |TaxonNamespace| and |Taxon| objects, which is important if you are going to be carrying comparisons or operations between multiple |Tree| objects.
As with the "|get_from_methods|" methods, keyword arguments can be used to provide :ref:`control on the data source parsing <Customizing_Data_Creation_and_Reading>`.


.. note:: DendroPy 3.xx supported "|read_from_methods|" methods on |Tree| and |CharacterMatrix|-derived classes. This is no longer supported in DendroPy 4 and above. Instead of trying to re-populate an existing |Tree| or |CharacterMatrix|-derived object by using "|read_from_methods|"::

            x = dendropy.Tree()
            x.read_from_path("tree1.nex", "nexus")
            .
            .
            .
            x.read_from_path("tree2.nex", "nexus")

        simply rebind the new object returned by "|get_from_methods|"::

            x = dendropy.Tree.get_from_path("tree1.nex", "nexus")
            .
            .
            .
            x = dendropy.Tree.get_from_path("tree2.nex", "nexus")

.. _Specifying_the_Data_Source_Format:

Specifying the Data Source Format
==================================

All the "|get_from_methods|" and "|read_from_methods|" methods take a :term:`schema` specification string using the ``schema`` argument which specifies the format of the data source.

The string can be one of the following:

    "``nexus``"
        To read |Tree|, |TreeList|, |CharacterMatrix|, or |DataSet| objects from a NEXUS-formatted source.

    "``newick``"
        To read |Tree|, |TreeList|, or |DataSet| objects from a Newick-formatted source.

    "``fasta``"
        To read |CharacterMatrix| or |DataSet| objects from a FASTA-formatted source. FASTA-sources require the additional keyword, ``data_type``, that describes the type of data: "``dna``", "``rna``", "``protein``", "``standard``"" (discrete data represented as binary 0/1), "``restriction``" (restriction sites), or "``infinite``" (infinite sites).

    "``phylip``"
        To read |CharacterMatrix| or |DataSet| objects from a PHYLIP-formatted source.
        You would typically use a specific |CharacterMatrix| class depending on the data type: e.g. |DnaCharacterMatrix|, |ContinuousCharacterMatrix| etc. If you use a more general class, e.g. |DataSet|, then for PHYLIP-sources you need to specify the additional keyword argument, ``data_type``, that describes the type of data: "``dna``", "``rna``", "``protein``", "``standard``"" (discrete data represented as binary 0/1), "``restriction``" (restriction sites), or "``infinite``" (infinite sites).

    "``beast-summary-tree``"
        To read |Tree| or |TreeList| objects from a BEAST annotated consensus tree source.
        Each node on the resulting tree(s) will have the following attributes: "``height``", "``height_median``", "``height_95hpd``", "``height_range``", "``length``", "``length_median``", "``length_95hpd``", "``length_range``", "``posterior'. Scalar values will be of ``float`` type, while ranges (e.g., "``height_95hpd``", "``height_range``", "``length_95hpd``", "``length_range``") will be two-element lists of ``float``.

.. _Customizing_Data_Creation_and_Reading:

Customizing Data Creation and Reading
=====================================

When specifying a data source from which to create or populate data objects
using the "|get_from_methods|" or "|read_from_methods|" methods, you can also
specify keyword arguments that provide fine-grained control over how the data
source is parsed.

Some of these keyword arguments apply generally, regardless of the format of
the data source or the data object being created, while others are specific to
the data object type or the data source format.

All Schemas
^^^^^^^^^^^

    ``attached_taxon_namespace``
        If |True| when reading into a |DataSet| object, then a new
        |TaxonNamespace| object will be created and added to the
        :attr:`~dendropy.datamodel.datasetmodel.DataSet.taxon_namespaces` list
        of the |DataSet| object, and the |DataSet| object will be placed in
        "attached" (or single) taxon set mode, i.e., all taxa in any data
        sources parsed or read will be mapped to the same |TaxonNamespace|
        object. By default, this is |False|, resulting in a multi-taxon set
        mode |DataSet| object.

    ``taxon_namespace``
        If passed a |TaxonNamespace| object, then this |TaxonNamespace| will be
        used to manage all taxon references in the data source.  When creating
        a new |Tree|, |TreeList| or |CharacterMatrix| object from a data
        source, the |TaxonNamespace| object passed by this keyword will be used
        as the |TaxonNamespace| associated with the object.
        When reading into a |DataSet| object, if the data source defines
        multiple collections of taxa (as is possible with, for example, the
        NEXML :term:`schema`, or the Mesquite variant of the NEXUS :term:`schema`), then
        multiple new |TaxonNamespace| object will be created. By passing a
        |TaxonNamespace| object through the ``taxon_namespace`` keyword, you
        can force DendroPy to use the same |TaxonNamespace| object for all
        taxon references.

    ``exclude_trees``
        If |True|, then all tree data in the data source will be skipped.
        Default value is |False|, i.e., all tree data will be included.

    ``exclude_chars``
        If |True|, then all character data in the data source will be skipped.
        Default value is |False|, i.e., all character data will be included.

Schema-Specific
^^^^^^^^^^^^^^^


Newick
......

.. autodocstringonly:: dendropy.dataio.newickreader.NewickReader.__init__

NEXUS
.....

.. autodocstringonly:: dendropy.dataio.nexusreader.NexusReader.__init__

FASTA
.....

.. autodocstringonly:: dendropy.dataio.fastareader.FastaReader.__init__

PHYLIP
......


BEAST Summary Trees
...................

