*************************************
Reading and Writing Phylogenetic Data
*************************************

Creating New Objects From an External Data Source
=================================================

The |Tree|, |TreeList|, |CharacterMatrix|-derived (i.e., |DnaCharacterMatrix|,
|ProteinCharacterMatrix|, |StandardCharacterMatrix|, etc.), and |DataSet|
classes all support a "|get|" factory class-method that instantiates an object
of the given class from a data source. This method takes, at a minumum, two
keyword arguments that specify the *source* of the data and the *schema* (or
format) of the data.

The source must be specifed using one and exactly one of the following:

    -   a path to a file (specified using the keyword argument "``path``")
    -   a file or a file-like object opened for reading (specified using the keyword argument ``"file"``)
    -   a string value giving the data directly (specified using the keyword argument ``"data"``)
    -   or a URL (specified using the keyword argument ``"url"``)

The schema is specified using the keyword argument ``"schema"``, and takes a string value that identifies the format of data.
This ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" can be one of: ":doc:`fasta </schemas/fasta>`", ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", ":doc:`nexml </schemas/nexml>`", or ":doc:`phylip </schemas/phylip>`".
Not all formats are supported for reading, and not all formats make sense for particular objects (for example, it would not make sense to try and instantiate a |Tree| or |TreeList| object from a |FASTA|-formatted data source).

.. A ":term:`schema`" is DendroPy-speak for "format" (we cannot use the argument or variable name "format" for this in library, because this is a Python built-in, and hence we use "schema" and adopted this terminology for consistency), and is specified using one of a set of predefined string values.

For example:

.. code-block:: python

    import dendropy

    tree1 = dendropy.Tree.get(path="mle.tre", schema="newick")
    tree2 = dendropy.Tree.get(file=open("mle.nex", "r"), schema="nexus")
    tree3 = dendropy.Tree.get(data="((A,B),(C,D));", schema="newick")
    tree4 = dendropy.Tree.get(url="http://api.opentreeoflife.org/v2/study/pg_1144/tree/tree2324.nex", schema="nexus")

    tree_list1 = dendropy.TreeList.get(path="pythonidae.mcmc.nex", schema="nexus")
    tree_list2 = dendropy.TreeList.get(file=open("pythonidae.mcmc.nex", "r"), schema="nexus")
    tree_list3 = dendropy.TreeList.get(data="(A,(B,C));((A,B),C);", "r"), schema="newick")

    dna1 = dendropy.DnaCharacterMatrix.get(file=open("pythonidae.fasta"), schema="fasta")
    dna2 = dendropy.DnaCharacterMatrix.get(url="http://purl.org/phylo/treebase/phylows/matrix/TB2:M2610?format=nexus", schema="nexus")
    aa1 = dendropy.ProteinCharacterMatrix.get(file=open("pythonidae.dat"), schema="phylip")
    std1 = dendropy.StandardCharacterMatrix.get(path="python_morph.nex", schema="nexus")
    std2 = dendropy.StandardCharacterMatrix.get(data=">t1\n01011\n\n>t2\n11100", schema="fasta")

    dataset1 = dendropy.DataSet.get(path="pythonidae.chars_and_trees.nex", schema="nexus")
    dataset2 = dendropy.DataSet.get(url="http://purl.org/phylo/treebase/phylows/study/TB2:S1925?format=nexml", schema="nexml")

The "|get|" method takes a number of other optional keyword arguments that provide control over how the data is interpreted and processed.
Some are general to all classes (e.g., the "``label``" or "``taxon_namespace``" arguments), while others specific to a given class (e.g. the "``exclude_trees``" argument when instantiating data into a |DataSet| object, or the "``tree_offset``" argument when instantiating data into a |Tree| or |TreeList| object).
These are all covered in detail in the documentation of the respective methods for each class:

    -   :meth:`Tree.get <dendropy.datamodel.treemodel.Tree.get>`
    -   :meth:`TreeList.get <dendropy.datamodel.treecollectionmodel.TreeList.get>`
    -   :meth:`DnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.DnaCharacterMatrix.get>`
    -   :meth:`RnaCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RnaCharacterMatrix.get>`
    -   :meth:`ProteinCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.ProteinCharacterMatrix.get>`
    -   :meth:`RestrictionSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterMatrix.get>`
    -   :meth:`InfiniteSitesCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterMatrix.get>`
    -   :meth:`StandardCharacterMatrix.get <dendropy.datamodel.charmatrixmodel.StandardCharacterMatrix.get>`
    -   :meth:`DataSet.get <dendropy.datamodel.datasetmodel.DataSet.get>`

Other optional keyword arguments are :ref:`specific to the schema or format <Schema_Specific_Keyword_Arguments>` (e.g., the "``preserve_underscores``" argument when reading |Newick| or |Nexus| data).
These are covered in detail in the :doc:`DendroPy Schema Guide </schemas/index>`.

.. note::

    The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes
    also support a "|get_from_methods|" family of factory class-methods that
    can be seen as specializations of the "|get|" method for various types of
    sources (in fact, the "|get|" method is actually a dispatcher that calls on
    one of these methods below for implementation of the functionality):

        :meth:`get_from_stream(src, schema, \*\*kwargs)`
            Takes a file or file-like object opened for reading the data source as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`get(file=src, schema=schema, ...)`".

        :meth:`get_from_path(src, schema, \*\*kwargs)`
            Takes a string specifying the path to the the data source file as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`get(path=src, schema=schema, ...)`".

        :meth:`get_from_string(src, schema, \*\*kwargs)`
            Takes a string containing the source data as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`get(data=src, schema=schema, ...)`".

        :meth:`get_from_url(src, schema, \*\*kwargs)`
            Takes a string containing the URL of the data as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.
            This is equivalent to calling ":meth:`get(url=src, schema=schema, ...)`".

    As with the "|get|" method, the additional keyword arguments are specific to the given class or schema type.

Adding Data to Existing Objects from an External Data Source
============================================================

In addition to the "|get|" class factory method, the collection classes (|TreeList|, |TreeArray| and |DataSet|) each support a "|read|" *instance* method that *add* data from external sources to an existing object (as opposed to creating and returning a new object based on an external data source).
This "|read|" instance method has a signature that parallels the "|get|" factory method described above, requiring:

    -   A specification of a source using one and exactly one of the following keyword arguments: "``path``", "``file``", "``data``", "``url``".
    -   A specification of the :ref:`schema <Specifying_the_Data_Source_Format>` or format of the data.
    -   Optional keyword arguments to customize/control the parsing and interpretation of the data.

As with the "|get|" method, the "|read|" method takes a number of other optional keyword arguments that provide control over how the data is interpreted and processed, which are covered in more detail in the documentation of the respective methods for each class:

    -   :meth:`TreeList.read <dendropy.datamodel.treecollectionmodel.TreeList.read>`
    -   :meth:`TreeArray.read <dendropy.datamodel.treecollectionmodel.TreeArray.read>`
    -   :meth:`DataSet.read <dendropy.datamodel.datasetmodel.DataSet.read>`

as well as :ref:`schema-specific keyword arguments <Schema_Specific_Keyword_Arguments>` which are covered in detail in the :doc:`DendroPy Schema Guide </schemas/index>`.

For example, the following accumulates post-burn-in trees from several different files into a single |TreeList| object::

    >>> import dendropy
    >>> post_trees = dendropy.TreeList()
    >>> post_trees.read(
    ...         file=open("pythonidae.nex.run1.t", "r")
    ...         schema="nexus",
    ...         tree_offset=200)
    >>> print(len(post_trees))
    800
    >>> post_trees.read(
    ...         path="pythonidae.nex.run2.t",
    ...         schema="nexus",
    ...         tree_offset=200)
    >>> print(len(post_trees))
    1600
    >>> s = open("pythonidae.nex.run3.t", "r").read()
    >>> post_trees.read(
    ...         data=s,
    ...         schema="nexus",
    ...         tree_offset=200)
    >>> print(len(post_trees))
    2400

.. The |TreeList| object automatically handles taxon management, and ensures that all appended |Tree| objects share the same |TaxonNamespace| reference. Thus all the |Tree| objects created and aggregated from the data sources in the example will all share the same |TaxonNamespace| and |Taxon| objects, which is important if you are going to be carrying comparisons or operations between multiple |Tree| objects.
.. As with the "|get|" method, keyword arguments can be used to provide :ref:`control on the data source parsing <Customizing_Data_Creation_and_Reading>`.

while the following accumulates data from a variety of sources into a single |DataSet| object under the same |TaxonNamespace| to ensure that they all reference the same set of |Taxon| objects::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> tns = ds.new_taxon_namespace()
    >>> ds.attach_taxon_namespace(tns)
    >>> ds.read(url="http://api.opentreeoflife.org/v2/study/pg_1144/tree/tree2324.nex",
    ...     schema="nexus")
    >>> ds.read(file=open("pythonidae.fasta"), schema="fasta")
    >>> ds.read(url="http://purl.org/phylo/treebase/phylows/matrix/TB2:M2610?format=nexus",
    ...     schema="nexus")
    >>> ds.read(file=open("pythonidae.dat"), schema="phylip")
    >>> ds.read(path="python_morph.nex", schema="nexus")
    >>> ds.read(data=">t1\n01011\n\n>t2\n11100", schema="fasta")

.. note:: DendroPy 3.xx supported "|read_from_methods|" methods on |Tree| and |CharacterMatrix|-derived classes. This is no longer supported in DendroPy 4 and above. Instead of trying to re-populate an existing |Tree| or |CharacterMatrix|-derived object by using "|read_from_methods|"::

            x = dendropy.Tree()
            x.read_from_path("tree1.nex", "nexus")
            .
            .
            .
            x.read_from_path("tree2.nex", "nexus")

        simply rebind the new object returned by "|get|"::

            x = dendropy.Tree.get(path="tree1.nex", schema="nexus")
            .
            .
            .
            x = dendropy.Tree.get(path="tree2.nex", schema="nexus")

.. note:: The |TreeList|, |TreeArray|, and |DataSet| classes
    also support a "|read_from_methods|" family of instance methods that
    can be seen as specializations of the "|read|" method for various types of
    sources (in fact, the "|read|" method is actually a dispatcher that calls on
    one of these methods below for implementation of the functionality):

        :meth:`read_from_stream(src, schema, \*\*kwargs)`
            Takes a file or file-like object opened for reading the data source as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`read(file=src, schema=schema, ...)`".

        :meth:`read_from_path(src, schema, \*\*kwargs)`
            Takes a string specifying the path to the the data source file as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`read(path=src, schema=schema, ...)`".

        :meth:`read_from_string(src, schema, \*\*kwargs)`
            Takes a string containing the source data as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be to control the parsing and other options.
            This is equivalent to calling ":meth:`read(data=src, schema=schema, ...)`".

        :meth:`read_from_url(src, schema, \*\*kwargs)`
            Takes a string containing the URL of the data as the first argument, and a :ref:`schema specification string <Specifying_the_Data_Source_Format>` as the second.
            Optional :term:`schema`-specific keyword arguments can be  to control the parsing and other options.
            This is equivalent to calling ":meth:`read(url=src, schema=schema, ...)`".

    As with the "|read|" method, the additional keyword arguments are specific to the given class or schema type.


Writing Out Phylogenetic Data
=============================

The |Tree|, |TreeList|, |CharacterMatrix|-derived (i.e., |DnaCharacterMatrix|,
|ProteinCharacterMatrix|, |StandardCharacterMatrix|, etc.), and |DataSet|
classes all support a "|write|" instance method for serialization of data to an
external data source.
This method takes two mandatory keyword arguments:

    -   One and exactly one of the following to specify the *destination*:
        -   a path to a file (specified using the keyword argument "``path``")
        -   a file or a file-like object opened for writing (specified using the keyword argument ``"file"``)

    -   A ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" given by the keyword argument "``schema``", to identify the schema or format for the output.

Alternatively, the |Tree|, |TreeList|, |CharacterMatrix|-derived, or |DnaCharacterMatrix| objects may also be represented as a string by calling the "``as_string()``" method, which requires at least one single mandatory argument, "``schema``", giving the ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" to identify the format of the output.

In either case, the ":ref:`schema specification string <Specifying_the_Data_Source_Format>`" can be one of: ":doc:`fasta </schemas/fasta>`", ":doc:`newick </schemas/newick>`", ":doc:`nexus </schemas/nexus>`", ":doc:`nexml </schemas/nexml>`", or ":doc:`phylip </schemas/phylip>`".

For example:

.. code-block:: python


    tree.write(path="output.tre", schema="newick")
    dest = open("output.xml", "w")
    tree_list.write(file=dest, schema="nexml")
    print(dna_character_matrix.as_string(schema="fasta"))


As with the "|get|" and "|read|" methods, further keyword arguments can be specified to control behavior.
These are covered in detail in the ":doc:`/schemas/index`" section.

.. note::

    The |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes also support a "|write_to_methods|" family of instance methods that can be seen as specializations of the "|write|" method for various types of destinations:

        :meth:`write_to_stream(dest, schema, \*\*kwargs)`
            Takes a file or file-like object opened for writing the data as the first argument, and a string specifying the schema as the second.

        :meth:`write_to_path(dest, schema, \*\*kwargs)`
            Takes a string specifying the path to the file as the first argument, and a string specifying the schema as the second.

        :meth:`as_string(schema, \*\*kwargs)`
            Takes a string specifying the schema as the first argument, and returns a string containing the formatted-representation of the data.

