*************************
Reading Phylogenetic Data
*************************

Creating and Populating New Objects
===================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes all support ":meth:`get_from_*()`" factory methods that allow for the simultaneous instantiation and population of the objects from a data source:

    :meth:`get_from_stream(src, schema, **kwargs)`
        Takes a file or file-like object opened for reading the data source as the first argument, and a string specifying the schema as the second.

    :meth:`get_from_path(src, schema, **kwargs)`
        Takes a string specifying the path to the the data source file as the first argument, and a string specifying the schema as the second.

    :meth:`get_from_string(src, schema, **kwargs)`
        Takes a string containing the source data as the first argument, and a string specifying the schema as the second.

A "schema" is DendroPy-speak for "format" (we cannot use the argument name "format" because this is a Python built-in, and hence we adopted this terminology for consistency), and is specified using one of a set of predefined string values.
The :ref:`schema specification string <Specifying_the_Data_Source_Format>` can be one of: "``nexus``", "``newick``", "``nexml``", "``fasta``", or "``phylip``".
Not all formats are supported for reading, and not all formats make sense for particular objects (for example, it would not make sense to try and instantiate a |Tree| or |TreeList| object from a FASTA-formatted data source).

All ":meth:`get_from_*()`"  methods minimally take a source and :ref:`schema specification string <Specifying_the_Data_Source_Format>` as arguments and return a new object of the given type populated from the given source::

    >>> import dendropy
    >>> tree1 = dendropy.Tree.get_from_string("((A,B),(C,D))", schema="newick")
    >>> tree_list1 = dendropy.TreeList.get_from_path("pythonidae.mcmc.nex", schema="nexus")
    >>> dna1 = dendropy.DnaCharacterMatrix.get_from_stream(open("pythonidae.fasta"), "dnafasta")
    >>> std1 = dendropy.StandardCharacterMatrix.get_from_path("python_morph.nex", "nexus")
    >>> dataset1 = dendropy.DataSet.get_from_path("pythonidae.nex", "nexus")

Various :ref:`keyword arguments <Customizing_Data_Creation_and_Reading>` can also be passed to these methods which customize or control how the data is parsed and mapped into DendroPy object space. These are discussed :ref:`below <Customizing_Data_Creation_and_Reading>`.

Reading and Populating (or Repopulating) Existing Objects
=========================================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes all support a suite of ":meth:`read_from_*()`" instance methods that parallels the ":meth:`get_from_*()`" factory methods described above:

    :meth:`read_from_stream(src, schema, **kwargs)`
        Takes a file or file-like object opened for reading the data source as the first argument, and a string specifying the schema as the second.

    :meth:`read_from_path(src, schema, **kwargs)`
        Takes a string specifying the path to the the data source file as the first argument, and a string specifying the schema as the second.

    :meth:`read_from_string(src, schema, **kwargs)`
        Takes a string specifying containing the source data as the first argument, and a string specifying the schema as the second.

When called on an existing |TreeList| or |DataSet| object, these methods *add* the data from the data source to the object, whereas when called on an existing |Tree| or |CharacterMatrix| object,  they *replace* the object's data with data from the data source.
As with the ":meth:`get_from_*()`" methods, the :ref:`schema specification string <Specifying_the_Data_Source_Format>` can be any supported and type-apppropriate schema, such as "``nexus``", "``newick``", "``nexml``", "``fasta``", "``phylip``", etc.

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

In contrast to the aggregating behavior of :meth:`read_from_*()` of |TreeList| and |DataSet| objects, the :meth:`read_from_*()` methods of |Tree|- and |CharacterMatrix|-derived objects show replacement behavior. For example, the following changes the contents of a |Tree| by re-reading it::

    >>> import dendropy
    >>> t = dendropy.Tree()
    >>> t.read_from_path('pythonidae.mle.nex', 'nexus')
    >>> print(t.description())
    Tree object at 0x79c70 (Tree37413776: '0'): ('Python molurus':0.0779719244,(('Python sebae':0.1414715009,((((('Morelia tracyae':0.0435011998,('Morelia amethistina':0.0305993564,(('Morelia nauta':0.0092774432,'Morelia kinghorni':0.0093145395):0.005595,'Morelia clastolepis':0.005204698):0.023435):0.012223):0.025359,'Morelia boeleni':0.0863199106):0.019894,(('Python reticulatus':0.0828549023,'Python timoriensis':0.0963051344):0.072003,'Morelia oenpelliensis':0.0820543043):0.002785):0.00274,(((('Morelia viridis':0.0925974416,('Morelia carinata':0.0943697342,('Morelia spilota':0.0237557178,'Morelia bredli':0.0357358071):0.041377):0.005225):0.004424,('Antaresia maculosa':0.1141193265,(('Antaresia childreni':0.0363195704,'Antaresia stimsoni':0.0188535952):0.043287,'Antaresia perthensis':0.0947695442):0.019148):0.007921):0.022413,('Leiopython albertisii':0.0698883547,'Bothrochilus boa':0.0811607602):0.020941):0.007439,(('Liasis olivaceus':0.0449896545,('Liasis mackloti':0.0331564496,'Liasis fuscus':0.0230286886):0.058253):0.016766,'Apodora papuana':0.0847328612):0.008417):0.006539):0.011557,('Aspidites ramsayi':0.0349772256,'Aspidites melanocephalus':0.0577536309):0.042499):0.036177):0.016859,'Python brongersmai':0.1147218285):0.001271,'Python regius':0.1800489093):0.0
    >>> t.read_from_path('pythonidae.mcmc-con.nex', 'nexus')
    >>> print(t.description())
    Tree object at 0x79c70 (Tree37414064: 'con 50 majrule'): ('Python regius':0.212275,('Python sebae':0.176816,(((((('Antaresia maculosa':0.127351,('Antaresia perthensis':0.108378,('Antaresia stimsoni':0.021372,'Antaresia childreni':0.038155):0.046446):0.025262):0.012957,('Morelia carinata':0.101145,('Morelia bredli':0.038563,'Morelia spilota':0.025643):0.050967):0.010472,'Morelia viridis':0.098541):0.023291,('Bothrochilus boa':0.091928,'Leiopython albertisii':0.080986):0.031583):0.008347,((('Liasis fuscus':0.026601,'Liasis mackloti':0.034524):0.069881,'Liasis olivaceus':0.047727):0.023758,'Apodora papuana':0.096097):0.01474):0.010084,(('Python timoriensis':0.101865,'Python reticulatus':0.095018):0.0922,('Morelia boeleni':0.093309,('Morelia tracyae':0.04727,('Morelia amethistina':0.034936,(('Morelia nauta':0.011,'Morelia kinghorni':0.011198):0.006932,'Morelia clastolepis':0.008103):0.025987):0.017415):0.033886):0.027519,'Morelia oenpelliensis':0.092143):0.006779):0.018238,('Aspidites ramsayi':0.030898,'Aspidites melanocephalus':0.068553):0.049525):0.050607):0.023304,('Python brongersmai':0.132193,'Python molurus':0.08872):0.011466)

As with the :meth:`get_from_*()` methods, keyword arguments can be used to provide :ref:`control on the data source parsing <Customizing_Data_Creation_and_Reading>`.

.. _Specifying_the_Data_Source_Format:

Specifying the Data Source Format
==================================

All the :meth:`get_from_*()` and :meth:`read_from_*()` methods take a schema specification string using the ``schema`` argument which specifies the format of the data source.

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
using the :meth:`get_from_*()`, :meth:`read_from_*()`, or passing a data source
stream to a constructor, you can also specify keyword arguments that provide
fine-grained control over how the data source is parsed.

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
        NEXML schema, or the Mesquite variant of the NEXUS schema), then
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

