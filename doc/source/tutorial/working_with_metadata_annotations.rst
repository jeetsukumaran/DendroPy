*************************************
Working with Metadata and Annotations
*************************************

|DendroPy| provides a rich infrastructure for decorating most types of phylogenetic objects (e.g., the |DataSet|, |TaxonSet|, |Taxon| |TreeList|, |Tree|, and various |CharacterMatrix| classes) with metadata information.
These phylogenetic objects have an attribute, :attr:`annotations`, that is an instance of the :class:`~dendropy.dataobject.base.AnnotationSet` class, which is an iterable (derived from :class:`set`) that serves to manage a collection of :class:`~dendropy.dataobject.base.Annotation` objects.
Each :class:`~dendropy.dataobject.base.Annotation` object tracks a single annotation element.
These annotations will be rendered as ``meta`` elements when writing to NeXML format or "hot comments" when writing to NEXUS/NEWICK format.
Note that full and robust expression of metadata annotations, including stable and consistent round-tripping of information, can only be achieved while in the NeXML format.

Basic Metadata Annotation Creation
==================================

The :meth:`~dendropy.dataobject.base.AnnotationSet.add_new` method of the :attr:`annotations` attribute allows for direct adding of metadata. This method has two mandatory arguments, ``name`` and ``value``::

    >>> import dendropy
    >>> tree = dendropy.Tree.get_from_path('pythonidae.mle.tree', 'nexus')
    >>> tree = dendropy.Tree.get_from_path('examples/pythonidae.mle.nex', 'nexus')
    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... )

When printing the tree in NeXML, the metadata will be rendered as a ``<meta>`` tag child element of the associated ``<tree>`` element::

    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009"
        xmlns:dendropy="http://packages.python.org/DendroPy/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
    >
    .
    .
    .
        <trees id="x4320340992" otus="x4320340552">
            <tree id="x4320381904" label="0" xsi:type="nex:FloatTree">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:subject" content="Python phylogenetics" id="meta4320379536" />
    .
    .
    .

As can be seen, by default the metadata property is mapped to the "``dendropy``" namespace (i.e., '``xmlns:dendropy="http://packages.python.org/DendroPy/"``').

This can be customized by using the "``name_prefix``" and "``namespace``" arguments to the call to :meth:`~dendropy.dataobject.base.AnnotationSet.add_new`::

    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... name_prefix="dc",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... )

This will result in the following NeXML fragment::

    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:dendropy="http://packages.python.org/DendroPy/"
    >
    .
    .
    .
        <trees id="x4320340904" otus="x4320340464">
            <tree id="x4320377872" label="0" xsi:type="nex:FloatTree">
                <meta xsi:type="nex:LiteralMeta" property="dc:subject" content="Python phylogenetics" id="meta4320375440" />
    .
    .
    .

Note that if the "``name_prefix``" or "``namespace``" must be specified simultaneously; that is, if one is specified, then the other must be specified as well.

For NeXML output, you can also specify a datatype::


    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... datatype_hint="xsd:string",
    ... )
    >>> tree.annotations.add_new(
    ... name="answer",
    ... value=42,
    ... datatype_hint="xsd:integer",
    ... )

When writing to NeXML, this will result in the following fragment::

    <trees id="x4320340992" otus="x4320340552">
        <tree id="x4320381968" label="0" xsi:type="nex:FloatTree">
            <meta xsi:type="nex:LiteralMeta" property="dendropy:answer" content="42" datatype="xsd:integer" id="meta4320379536" />
            <meta xsi:type="nex:LiteralMeta" property="dendropy:subject" content="Python phylogenetics" datatype="xsd:string" id="meta4320379472" />

You can also specify that the data should be interpreted as a source to be dereferenced in NeXML by passing in ``compose_as_reference=True``.
Note that this does not actually populated the contents of the annotation from the source (unlike the dynamic attribute value binding discussed below), but just indicates the the contents of the annotation should be *interpreted* differently.
Thus, the following annotation::

    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="http://en.wikipedia.org/wiki/Pythonidae",
    ... name_prefix="dc",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... compose_as_reference=True,
    ... )

will be rendered in NeXML as::

    <meta xsi:type="nex:ResourceMeta" rel="dc:subject" href="http://en.wikipedia.org/wiki/Pythonidae" />
