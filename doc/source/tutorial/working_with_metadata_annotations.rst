*********************************
Working with Metadata Annotations
*********************************

|DendroPy| provides a rich infrastructure for decorating most types of phylogenetic objects (e.g., the |DataSet|, |TaxonNamespace|, |Taxon| |TreeList|, |Tree|, and various |CharacterMatrix| classes) with metadata information.
These phylogenetic objects have an attribute, :attr:`annotations`, that is an instance of the :class:`~dendropy.datamodel.basemodel.AnnotationSet` class, which is an iterable (derived from :class:`dendropy.utility.containers.OrderedSet`) that serves to manage a collection of :class:`~dendropy.datamodel.basemodel.Annotation` objects.
Each :class:`~dendropy.datamodel.basemodel.Annotation` object tracks a single annotation element.
These annotations will be rendered as ``meta`` elements when writing to NeXML format or ampersand-prepended comemnt strings when writing to NEXUS/NEWICK format.
Note that full and robust expression of metadata annotations, including stable and consistent round-tripping of information, can only be achieved while in the NeXML format.

Overview of the Infrastructure for Metadata Annotation in |DendroPy|
====================================================================

Each item of metadata is maintained in an object of the :class:`~dendropy.datamodel.basemodel.Annotation` class.
This class has the following attributes:

    :attr:`~dendropy.datamodel.basemodel.Annotation.name`
        The name of the metadata item or annotation.

    :attr:`~dendropy.datamodel.basemodel.Annotation.value`
        The value or content of the metadata item or annotation.

    :attr:`~dendropy.datamodel.basemodel.Annotation.datatype_hint`
        Custom data type indication for NeXML output (e.g. "xsd:string").

    :attr:`~dendropy.datamodel.basemodel.Annotation.name_prefix`
        Prefix that represents an abbreviation of the namespace associated with
        this metadata item.

    :attr:`~dendropy.datamodel.basemodel.Annotation.namespace`
        The namespace (e.g. "http://www.w3.org/XML/1998/namespace") of this
        metadata item (NeXML output).

    :attr:`~dendropy.datamodel.basemodel.Annotation.annotate_as_reference`
        If |True|, indicates that this annotation should not be interpreted semantically as a literal value, but rather as a source to be dereferenced.

    :attr:`~dendropy.datamodel.basemodel.Annotation.is_hidden`
        If |True|, indicates that this annotation should not be printed or written out.

    :attr:`~dendropy.datamodel.basemodel.Annotation.prefixed_name`
        Returns the name of this annotation with its namespace prefix (e.g. "dc:subject").

These :class:`~dendropy.datamodel.basemodel.Annotation` objects are typically collected and managed in a "annotations manager" container class, :class:`~dendropy.datamodel.basemodel.AnnotationSet`.
This is a specialization of :class:`dendropy.utility.containers.OrderedSet` whose elements are instances of :class:`~dendropy.datamodel.basemodel.Annotation`.
The full set of annotations associated with each object of |DataSet|, |TaxonNamespace|, |Taxon| |TreeList|, |Tree|, various |CharacterMatrix| and other phylogenetic data class types is available through the :attr:`annotations` attribute of those objects, which is an instance of :class:`~dendropy.datamodel.basemodel.AnnotationSet`.
The :class:`~dendropy.datamodel.basemodel.AnnotationSet` includes the following additional methods to support the creation, access, and management of the :class:`~dendropy.datamodel.basemodel.Annotation` object elements contained within it:

    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_new()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_bound_attribute()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_citation()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.find()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.drop()`
    - :meth:`~dendropy.datamodel.basemodel.AnnotationSet.values_as_dict()`


..
    The fundamental unit of metadata in |DendroPy| is the :class:`~dendropy.datamodel.basemodel.Annotation` object.
    Each :class:`~dendropy.datamodel.basemodel.Annotation` object stores information regarding a single item of metadata, keeping track of, at a minimum, the name and value or content of the metadata item, which is accessible through the attributes ":attr:`~dendropy.datamodel.basemodel.Annotation.name`" and  ":attr:`~dendropy.datamodel.basemodel.Annotation.value`" respectively.
    These :class:`~dendropy.datamodel.basemodel.Annotation` objects are typically collected and managed in a "annotations manager" container class, :class:`~dendropy.datamodel.basemodel.AnnotationSet`, which is a specialization of ":class:`dendropy.utility.containers.OrderedSet`".
    Phylogenetic data objects of |DataSet|, |TaxonNamespace|, |Taxon| |TreeList|, |Tree|, various |CharacterMatrix| and other classes all have an attribute, ":attr:`annotations`", that represents an instance of the :class:`~dendropy.datamodel.basemodel.AnnotationSet` class, and whose elements are :class:`~dendropy.datamodel.basemodel.Annotation` objects that collectively make up the full set of annotations or metadata associated with that particular phylogenetic data object.
        The elements of the ":attr:`annotations`" attribute of phylogenetic data objects are objects of :class:`~dendropy.datamodel.basemodel.Annotation` that collectively make up the full set of annotations or metadata associated with that particular phylogenetic data object.


The following code snippet reads in a data file in NeXML format, and dumps out the annotations::

    #! /usr/bin/env python

    import sys
    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml", "nexml")
    print "-- (dataset) ---\n"
    for a in ds.annotations:
        print "%s = '%s'" % (a.name, a.value)
    for tree_list in ds.tree_lists:
        for tree in tree_list:
            print "\n-- (tree '%s') --\n" % tree.label
            for a in tree.annotations:
                print "%s = '%s'" % (a.name, a.value)

Running the above results in::

    -- (dataset) ---

    bibliographicCitation = 'Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G. 2012. Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific. Systematics and Biodiversity, .'
    subject = 'whale-fall'
    changeNote = 'Generated on Wed Jun 06 11:02:45 EDT 2012'
    subject = 'wood-fall'
    title = 'Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific'
    publicationName = 'Systematics and Biodiversity'
    creator = 'Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G.'
    publisher = 'Systematics and Biodiversity'
    contributor = 'Wiklund H.'
    volume = ''
    contributor = 'Altamira I.V.'
    number = ''
    contributor = 'Glover A.'
    historyNote = 'Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlDocumentWriter@645f9132 $Rev: 1060 $'
    contributor = 'Smith C.'
    modificationDate = '2012-06-04'
    contributor = 'Baco A.'
    contributor = 'Dahlgren T.G.'
    identifier.study.tb1 = 'None'
    publicationDate = '2012'
    section = 'Study'
    doi = ''
    title.study = 'Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific'
    subject = 'New species'
    subject = 'Ophryotrocha'
    creationDate = '2012-05-09'
    subject = 'polychaeta'
    date = '2012-06-04'
    subject = 'molecular phylogeny'
    identifier.study = '12713'

    -- (tree 'con 50 majrule') --

    ntax.tree = '41'
    kind.tree = 'Species Tree'
    quality.tree = 'Unrated'
    isDefinedBy = 'http://purl.org/phylo/treebase/phylows/study/TB2:S12713'
    type.tree = 'Consensus'


The following sections discuss these methods and attributes in detail, describing how the create, read, write, search, and manipulate annotations.

Metadata Annotation Creation
=============================

Reading Data from an External Source
------------------------------------

When reading data in NeXML format, metadata annotations given in the source are automatically created and associated with the corresponding data objects.

The metadata annotations associated with the phylogenetic data objects are collected in the attribute ``annotations`` of the objects, which is an object of type :class:`~dendropy.datamodel.basemodel.AnnotationSet`.
Each annotation item is represented as an
object of type :class:`~dendropy.datamodel.basemodel.Annotation`.

For example::

    #! /usr/bin/env python
    import dendropy
    ds = dendropy.DataSet.get_from_path("pythonidae.annotated.nexml",
    "nexml")
    for a in ds.annotations:
        print "Data Set '%s': %s" % (ds.label, a)
    for taxon_namespace in ds.taxon_namespaces:
        for a in taxon_namespace.annotations:
            print "Taxon Set '%s': %s" % (taxon_namespace.label, a)
        for taxon in taxon_namespace:
            for a in taxon.annotations:
                print "Taxon '%s': %s" % (taxon.label, a)
    for tree_list in ds.tree_lists:
        for a in tree_list.annotations:
            print "Tree List '%s': %s" % (tree_list.label, a)
        for tree in tree_list:
            for a in tree.annotations:
                print "Tree '%s': %s" % (tree.label, a)

produces::

    Data Set 'None': description="composite dataset of Pythonid sequences and trees"
    Data Set 'None': subject="Pythonidae"
    Taxon Set 'None': subject="Pythonidae"
    Taxon 'Python regius': closeMatch="http://purl.uniprot.org/taxonomy/51751"
    Taxon 'Python sebae': closeMatch="http://purl.uniprot.org/taxonomy/51752"
    Taxon 'Python molurus': closeMatch="http://purl.uniprot.org/taxonomy/51750"
    Taxon 'Python curtus': closeMatch="http://purl.uniprot.org/taxonomy/143436"
    Taxon 'Morelia bredli': closeMatch="http://purl.uniprot.org/taxonomy/461327"
    Taxon 'Morelia spilota': closeMatch="http://purl.uniprot.org/taxonomy/51896"
    Taxon 'Morelia tracyae': closeMatch="http://purl.uniprot.org/taxonomy/129332"
    Taxon 'Morelia clastolepis': closeMatch="http://purl.uniprot.org/taxonomy/129329"
    Taxon 'Morelia kinghorni': closeMatch="http://purl.uniprot.org/taxonomy/129330"
    Taxon 'Morelia nauta': closeMatch="http://purl.uniprot.org/taxonomy/129331"
    Taxon 'Morelia amethistina': closeMatch="http://purl.uniprot.org/taxonomy/51895"
    Taxon 'Morelia oenpelliensis': closeMatch="http://purl.uniprot.org/taxonomy/461329"
    Taxon 'Antaresia maculosa': closeMatch="http://purl.uniprot.org/taxonomy/51891"
    Taxon 'Antaresia perthensis': closeMatch="http://purl.uniprot.org/taxonomy/461324"
    Taxon 'Antaresia stimsoni': closeMatch="http://purl.uniprot.org/taxonomy/461325"
    Taxon 'Antaresia childreni': closeMatch="http://purl.uniprot.org/taxonomy/51888"
    Taxon 'Morelia carinata': closeMatch="http://purl.uniprot.org/taxonomy/461328"
    Taxon 'Morelia viridisN': closeMatch="http://purl.uniprot.org/taxonomy/129333"
    Taxon 'Morelia viridisS': closeMatch="http://purl.uniprot.org/taxonomy/129333"
    Taxon 'Apodora papuana': closeMatch="http://purl.uniprot.org/taxonomy/129310"
    Taxon 'Liasis olivaceus': closeMatch="http://purl.uniprot.org/taxonomy/283338"
    Taxon 'Liasis fuscus': closeMatch="http://purl.uniprot.org/taxonomy/129327"
    Taxon 'Liasis mackloti': closeMatch="http://purl.uniprot.org/taxonomy/51889"
    Taxon 'Antaresia melanocephalus': closeMatch="http://purl.uniprot.org/taxonomy/51883"
    Taxon 'Antaresia ramsayi': closeMatch="http://purl.uniprot.org/taxonomy/461326"
    Taxon 'Liasis albertisii': closeMatch="http://purl.uniprot.org/taxonomy/129326"
    Taxon 'Bothrochilus boa': closeMatch="http://purl.uniprot.org/taxonomy/461341"
    Taxon 'Morelia boeleni': closeMatch="http://purl.uniprot.org/taxonomy/129328"
    Taxon 'Python timoriensis': closeMatch="http://purl.uniprot.org/taxonomy/51753"
    Taxon 'Python reticulatus': closeMatch="http://purl.uniprot.org/taxonomy/37580"
    Taxon 'Xenopeltis unicolor': closeMatch="http://purl.uniprot.org/taxonomy/196253"
    Taxon 'Candoia aspera': closeMatch="http://purl.uniprot.org/taxonomy/51853"
    Taxon 'Loxocemus bicolor': closeMatch="http://purl.uniprot.org/taxonomy/39078"
    Tree '0': treeEstimator="RAxML"
    Tree '0': substitutionModel="GTR+G+I"

Metadata annotations in NEXUS and NEWICK must be given in the form of "hot comments" either in BEAST/FigTree syntax::

    [&subject='Pythonidae']

    [&length_hpd95={0.01917252,0.06241567},length_quant_5_95={0.02461821,0.06197141},length_range={0.01570374,0.07787249},length_mean=0.0418470252488,length_median=0.04091105,length_sd=0.0113086027131]

or NHX-like syntax::

    [&&subject='Pythonidae']

    [&&length_hpd95={0.01917252,0.06241567},length_quant_5_95={0.02461821,0.06197141},length_range={0.01570374,0.07787249},length_mean=0.0418470252488,length_median=0.04091105,length_sd=0.0113086027131]

However, by default these annotations are not parsed into |DendroPy| data model
unless the keyword argument ``extract_comment_metadata=True`` is passed in to the call::

    >>> ds = dendropy.DataSet.get_from_path("data.nex",
    ... "nexus",
    ... extract_comment_metadata=True)

In general, support for metadata in NEXUS and NEWICK formats is very basic and lossy, and is limited to a small range of phylogenetic data types (taxa, trees, nodes, edges).
These issues and limits are fundamental to the NEXUS and NEWICK formats, and thus if metadata is important to you and your work, you should be working with NeXML format.
The NeXML format provides for rich, flexible and robust metadata annotation for the broad range of phylogenetic data, and |DendroPy| provides full support for metadata reading and writing in NeXML.


Direct Composition with Literal Values
--------------------------------------

The :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_new` method of the :attr:`annotations` attribute allows for direct adding of metadata. This method has two mandatory arguments, "``name``" and "``value``"::

    >>> import dendropy
    >>> tree = dendropy.Tree.get_from_path('pythonidae.mle.tree', 'nexus')
    >>> tree = dendropy.Tree.get_from_path('examples/pythonidae.mle.nex', 'nexus')
    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... )

When printing the tree in NeXML, the metadata will be rendered as a "``<meta>``" tag child element of the associated "``<tree>``" element::

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

As can be seen, by default, the metadata property is mapped to the "``dendropy``" namespace (i.e., '``xmlns:dendropy="http://packages.python.org/DendroPy/"``').
This can be customized by using the "``name_prefix``" and "``namespace``" arguments to the call to :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_new`::

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

Note that the "``name_prefix``" or "``namespace``" must be specified simultaneously; that is, if one is specified, then the other must be specified as well.
For convenience, you can specify the name of the annotation with the name prefix prepended by specifying "``name_is_prefixed=True``", though the namespace must still be provided separately::

    >>> tree.annotations.add_new(
    ... name="dc:subject",
    ... value="Python phylogenetics",
    ... name_is_prefixed=True,
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... )

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

You can also specify that the data should be interpreted as a source to be dereferenced in NeXML by passing in ``annotate_as_reference=True``.
Note that this does not actually populate the contents of the annotation from the source (unlike the dynamic attribute value binding discussed below), but just indicates the the contents of the annotation should be *interpreted* differently by semantic readers.
Thus, the following annotation::

    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="http://en.wikipedia.org/wiki/Pythonidae",
    ... name_prefix="dc",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... annotate_as_reference=True,
    ... )

will be rendered in NeXML as::

    <meta xsi:type="nex:ResourceMeta" rel="dc:subject" href="http://en.wikipedia.org/wiki/Pythonidae" />

Sometimes, you may want to annotate an object with metadata, but do not want it to be printed or written out.
Passing the ``is_hidden=True`` argument will result in the annotation being suppressed in all output::

    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... name_prefix="dc",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... is_hidden=True,
    ... )

The ``is_hidden`` attribute of the an :class:`~dendropy.datamodel.basemodel.Annotation` object can also be set directly::

    >>> subject_annotations = tree.annotations.findall(name="citation")
    >>> for a in subject_annotations:
    ...    a.is_hidden = True

Dynamically Binding Annotation Values to Object Attribute Values
----------------------------------------------------------------

In some cases, instead of "hard-wiring" in metadata for an object, you may want to write out metadata that takes its value from the value of an attribute of the object.
The :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_bound_attribute` method allows you to do this.
This method takes, as a minimum, a *string* specifying the *name* of an existing attribute to which the value of the annotation will be dynamically bound.

For example:

.. literalinclude:: /examples/dynamic_annotations1.py

results in::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009"
        xmlns:dendropy="http://packages.python.org/DendroPy/"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
    >
        <otus id="x4320344648">
            <otu id="x4320380112" label="A">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:category" content="tiny" id="meta4320379472" />
            </otu>
            <otu id="x4320380432" label="B">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:category" content="medium" id="meta4320379536" />
            </otu>
            <otu id="x4320380752" label="C">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:category" content="N/A" id="meta4320379792" />
            </otu>
            <otu id="x4320381072" label="D">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:category" content="tiny" id="meta4320381328" />
            </otu>
            <otu id="x4320381264" label="E">
                <meta xsi:type="nex:LiteralMeta" property="dendropy:category" content="tiny" id="meta4320381392" />
            </otu>
        </otus>
        <trees id="x4320344560" otus="x4320344648">
            <tree id="x4320379600" xsi:type="nex:FloatTree">
                <node id="x4320379856">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="5491" id="meta4320379280" />
                </node>
                <node id="x4320379984" otu="x4320380112">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="2721" id="meta4320379408" />
                </node>
                <node id="x4320380176">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="4627" id="meta4320379344" />
                </node>
                <node id="x4320380304" otu="x4320380432">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="7202" id="meta4320381456" />
                </node>
                <node id="x4320380496">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="5337" id="meta4320379664" />
                </node>
                <node id="x4320380624" otu="x4320380752">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="1478" id="meta4320381520" />
                </node>
                <node id="x4320380816">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="1539" id="meta4320379728" />
                </node>
                <node id="x4320380944" otu="x4320381072">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="3457" id="meta4320381584" />
                </node>
                <node id="x4320381136" otu="x4320381264">
                    <meta xsi:type="nex:LiteralMeta" property="dendropy:pop_size" content="3895" id="meta4320381648" />
                </node>
                <rootedge id="x4320379920" target="x4320379856" />
                <edge id="x4320380048" source="x4320379856" target="x4320379984" />
                <edge id="x4320380240" source="x4320379856" target="x4320380176" />
                <edge id="x4320380368" source="x4320380176" target="x4320380304" />
                <edge id="x4320380560" source="x4320380176" target="x4320380496" />
                <edge id="x4320380688" source="x4320380496" target="x4320380624" />
                <edge id="x4320380880" source="x4320380496" target="x4320380816" />
                <edge id="x4320381008" source="x4320380816" target="x4320380944" />
                <edge id="x4320381200" source="x4320380816" target="x4320381136" />
            </tree>
        </trees>
    </nex:nexml>


By default, the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_bound_attribute` method uses the name of the attribute as the name of the annotation.
The "``annotation_name``" argument allows you explictly set the name of the annotation.
In addition, the method call also supports the other customization arguments of the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_new` method: "``datatype_hint``", "``name_prefix``", "``namespace``", "``name_is_prefixed``", "``annotate_as_reference``", "``is_hidden``", etc.::

    >>> tree.source_uri = None
    >>> tree.annotations.add_bound_attribute(
    ... "source_uri",
    ... annotation_name="dc:subject",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... annotate_as_reference=True)

Adding Citation Metadata
------------------------

You can add citation annotations using the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_citation` method.
This method takes at least one argument, ``citation``.
This can be a string representing the citation as a BibTex record or a dictionary with BibTex fields as keys and field content as values.

For example:

.. literalinclude:: /examples/bibtex_annotations1.py

will result in::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009"
        xmlns:bibtex="http://www.edutella.org/bibtex#"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:dendropy="http://packages.python.org/DendroPy/"
    >
        <meta xsi:type="nex:LiteralMeta" property="bibtex:journal" content="Molecular Biology and Evolution" datatype="xsd:string" id="meta4320453648" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:bibtype" content="article" datatype="xsd:string" id="meta4320453200" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:number" content="3" datatype="xsd:string" id="meta4320453776" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:citekey" content="HeathHH2012" datatype="xsd:string" id="meta4320453328" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:pages" content="939-955" datatype="xsd:string" id="meta4320453968" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:volume" content="29" datatype="xsd:string" id="meta4320453840" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:year" content="2012" datatype="xsd:string" id="meta4320453904" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:doi" content="10.1093/molbev/msr255" datatype="xsd:string" id="meta4320453456" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:title" content="A {Dirichlet} Process Prior for Estimating Lineage-Specific Substitution Rates." datatype="xsd:string" id="meta4320453520" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:url" content="http://mbe.oxfordjournals.org/content/early/2011/11/04/molbev.msr255.abstract" datatype="xsd:string" id="meta4320453584" />
        <meta xsi:type="nex:LiteralMeta" property="bibtex:author" content="Tracy A. Heath and Mark T. Holder and John P. Huelsenbeck" datatype="xsd:string" id="meta4320453712" />
        .
        .
        .


The following results in the same output as above, but the citation is given as a dictionary with BibTex fields as keys and content as values:

.. literalinclude:: /examples/bibtex_annotations2.py

By default, the citation gets annotated as a series of separate BibTex elements.
You can specify alternate formats by using the "``store_as``" argument.
This argument can take one of the following values:

    - "``bibtex``"
        Each BibTex field gets recorded as a separate annotation, with name
        given by the field name, content by the field value.
        This is the default, and the results in NeXML are shown above.

    - "``dublin``"
        A subset of the BibTex fields gets recorded as a set of Dublin Core (Publishing Requirements for Industry Standard Metadata) annotations, one per field::

        <meta xsi:type="nex:LiteralMeta" property="dc:date" content="2012" datatype="xsd:string" id="meta4320461584" />
        <meta xsi:type="nex:LiteralMeta" property="dc:publisher" content="Molecular Biology and Evolution" datatype="xsd:string" id="meta4320461648" />
        <meta xsi:type="nex:LiteralMeta" property="dc:title" content="A {Dirichlet} Process Prior for Estimating Lineage-Specific Substitution Rates." datatype="xsd:string" id="meta4320461776" />
        <meta xsi:type="nex:LiteralMeta" property="dc:creator" content="Tracy A. Heath and Mark T. Holder and John P. Huelsenbeck" datatype="xsd:string" id="meta4320461712" />

    - "``prism``"
        A subset of the BibTex fields gets recorded as a set of PRISM (Publishing Requirements for Industry Standard Metadata) annotations, one per field::

        <meta xsi:type="nex:LiteralMeta" property="prism:volume" content="29" datatype="xsd:string" id="meta4320461584" />
        <meta xsi:type="nex:LiteralMeta" property="prism:pageRange" content="939-955" datatype="xsd:string" id="meta4320461648" />
        <meta xsi:type="nex:LiteralMeta" property="prism:publicationDate" content="2012" datatype="xsd:string" id="meta4320461776" />
        <meta xsi:type="nex:LiteralMeta" property="prism:publicationName" content="Molecular Biology and Evolution" datatype="xsd:string" id="meta4320461712" />


In addition, the method call also supports some of the other customization arguments of the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.add_new` method:  "``name_prefix``", "``namespace``", "``name_is_prefixed``", "``is_hidden``".

Copying Metadata Annotations from One Phylogenetic Data Object to Another
-------------------------------------------------------------------------

As the :class:`~dendropy.datamodel.basemodel.AnnotationSet` is derived from :class:`dendropy.utility.containers.OrderedSet`, it has the :meth:`dendropy.utility.containers.OrderedSet.add` and :meth:`dendropy.utility.containers.OrderedSet.update` methods available for direct addition of :class:`~dendropy.datamodel.basemodel.Annotation` objects.
The following example shows how to add metadata annotations associated with a |DataSet| object to all its |Tree| objects::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    ds_annotes = ds.annotations.findall(name_prefix="dc").values_as_dict()
    for tree_list in ds.tree_lists:
        for tree in tree_list:
            tree.annotations.update(ds_annotes)

Or, alternatively::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    ds_annotes = ds.annotations.findall(name_prefix="dc").values_as_dict()
    for tree_list in ds.tree_lists:
        for tree in tree_list:
            for a in ds_annotes:
                tree.annotations.add(a)


Metadata Annotation Access and Manipulation
===========================================

Iterating Over Collections of Annotations
-----------------------------------------

The collection of :class:`~dendropy.datamodel.basemodel.Annotation` objects representing metadata annotations associated with particular phylgoenetic data objects can be accessed through the :attr:`annotations` attribute of each particular object.

For example::

    #! /usr/bin/env python
    ds = dendropy.DataSet.get_from_path("pythonidae.annotated.nexml",
    "nexml")
    for a in ds.annotations:
        print "The dataset has metadata annotation '%s' with content '%s'" % (a.name, a.value)
    tree = ds.tree_lists[0][0]
    for a in tree.annotations:
        print "Tree '%s' has metadata annotation '%s' with content '%s'" % (tree.label, a.name, a.value)

will result in::

    The dataset has metadata annotation 'description' with content 'composite dataset of Pythonid sequences and trees'
    The dataset has metadata annotation 'subject' with content 'Pythonidae'
    Tree '0' has metadata annotation 'treeEstimator' with content 'RAxML'
    Tree '0' has metadata annotation 'substitutionModel' with content 'GTR+G+I'

Retrieving Annotations By Search Criteria
-----------------------------------------

Instead of interating through every element in the :attr:`annotations` attribute of data objects, you can use the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` method of the the :attr:`annotations` object to return a *collection* of :class:`~dendropy.datamodel.basemodel.Annotation` objects that match the search or filter criteria specified in keyword arguments to the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` call.
These keyword arguments should specify attributes of :class:`~dendropy.datamodel.basemodel.Annotation` and the corresponding value to be matched.
Multiple keyword-value pairs can be specified, and only :class:`~dendropy.datamodel.basemodel.Annotation` objects that match *all* the criteria will be returned.

For example, the following returns a collection of annotations that have a name of "contributor"::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    results = ds.annotations.findall(name="contributor")
    for a in results:
        print "%s='%s'" % (a.name, a.value)

and will result in::

    contributor='Dahlgren T.G.'
    contributor='Baco A.'
    contributor='Smith C.'
    contributor='Glover A.'
    contributor='Altamira I.V.'
    contributor='Wiklund H.'

While the following returns a collection of annotations that are in the Dublin Core namespace::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    results = ds.annotations.findall(namespace="http://purl.org/dc/elements/1.1/")
    for a in results:
        print "%s='%s'" % (a.name, a.value)

and results in::

    subject='wood-fall'
    contributor='Wiklund H.'
    publisher='Systematics and Biodiversity'
    subject='whale-fall'
    contributor='Dahlgren T.G.'
    contributor='Smith C.'
    date='2012-06-04'
    subject='polychaeta'
    contributor='Glover A.'
    subject='Ophryotrocha'
    title='Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific'
    subject='New species'
    subject='molecular phylogeny'
    contributor='Altamira I.V.'
    creator='Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G.'
    contributor='Baco A.'

The following, in turn, searches for and suppresses printing of annotations that have a name prefix of "dc" *and* have empty values::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    results = ds.annotations.findall(name_prefix="dc", value="")
    for a in results:
        a.is_hidden = True

Modifying the :class:`~dendropy.datamodel.basemodel.Annotation` objects in a returned collection modifies the metadata of the parent data object. For example, the following sets all the field values to upper case characters::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    results = ds.annotations.findall(name="contributor")
    for a in results:
        a.value = a.value.upper()
    results = ds.annotations.findall(name="contributor")
    for a in results:
        print a.value

and results in::

    DAHLGREN T.G.
    BACO A.
    SMITH C.
    GLOVER A.
    ALTAMIRA I.V.
    WIKLUND H.

The collection returned by the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` method is an object of type :class:`~dendropy.datamodel.basemodel.AnnotationSet`.
However, while modifying :class:`~dendropy.datamodel.basemodel.Annotation` objects in this collection will result in the metadata of the parent object being modified (as in the previous example), adding new annotations to this returned collection will *not*  add them to the collection of metadata annotations of the parent object.
Thus, the following example shows that the size of the annotations collection associated with the dataset is unchanged by adding new annotations to the results of a :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` call::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    print len(ds.annotations)
    results = ds.annotations.findall(namespace="http://purl.org/dc/elements/1.1/")
    results.add_new(name="color", value="blue")
    results.add_new(name="height", value="100")
    results.add_new(name="length", value="200")
    results.add_new(name="width", value="50")
    print len(ds.annotations)

The above produces::

    30
    30

As can be seen, no new annotations are added to the data set metadata.

If no matching :class:`~dendropy.datamodel.basemodel.Annotation` objects are found then the :class:`~dendropy.datamodel.basemodel.AnnotationSet` that is returned is empty.

If *no* keyword arguments are passed to :meth:`~dendropy.datamodel.basemodel.Annotation.findall`, then *all* annotations are returned::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    results = ds.annotations.findall()
    print len(results) == len(ds.annotations)

The above produces::

    True

Retrieving a Single Annotation By Search Criteria
-------------------------------------------------

The :meth:`~dendropy.datamodel.basemodel.AnnotationSet.find` method of the the :attr:`annotations` object return a the *first* :class:`~dendropy.datamodel.basemodel.Annotation` object that matches the search or filter criteria specified in keyword arguments to the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` call.
These keyword arguments should specify attributes of :class:`~dendropy.datamodel.basemodel.Annotation` and the corresponding value to be matched.
Multiple keyword-value pairs can be specified, and only the first :class:`~dendropy.datamodel.basemodel.Annotation` object that matches *all* the criteria will be returned.

For example::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    print ds.annotations.find(name="contributor")

and will result in::

    contributor='Dahlgren T.G.'

While the following returns the first annotation in the Dublin Core namespace::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    print ds.annotations.find(namespace="http://purl.org/dc/elements/1.1/")

and results in::

    subject='wood-fall'

If no matching :class:`~dendropy.datamodel.basemodel.Annotation` objects are found then a default of |None| is returned::

    >>> print ds.annotations.find(name="author")
    None

Unlike :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall`, it is invalid to call :meth:`~dendropy.datamodel.basemodel.AnnotationSet.find` with no search criteria keyword arguments, and an ``TypeError`` exception will be raised.

Retrieving the Value of a Single Annotation
-------------------------------------------

For convenience, the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.get_value`, method is provided.
This will search the :class:`~dendropy.datamodel.basemodel.AnnotationSet` for the *first* :class:`~dendropy.datamodel.basemodel.Annotation` that has its name field equal to the first argument passed to the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.get_value` method, and return its value.
If no match is found, the second argument is returned (or |None|, if no second argument is specified).
Examples::

    >>> print tree.annotations.get_value("subject")
    molecular phylogeny
    >>> print tree.annotations.get_value("creator")
    Yoder A.D., & Yang Z.
    >>> print tree.annotations.get_value("generator")
    None
    >>> print tree.annotations.get_value("generator", "unspecified")
    unspecified

Transforming Annotations to a Dictionary
----------------------------------------

In some applications, it might be more convenient to work with dictionaries rather than :class:`~dendropy.datamodel.basemodel.AnnotationSet` objects.
The :meth:`~dendropy.datamodel.basemodel.Annotation.values_as_dict` methods creates a dictionary populated with key-value pairs from the collection.
By default, the keys are the ``name`` attribute of the :class:`~dendropy.datamodel.basemodel.Annotation` object and the values are the ``value`` attribute.
Thus, the following::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    a = ds.annotations.values_as_dict()
    print a

results in::

    {'volume': '',
    'doi': '',
    'date': '2012-06-04',
    'bibliographicCitation': 'Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G. 2012. Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific. Systematics and Biodiversity, .',
    'changeNote': 'Generated on Wed Jun 06 11:02:45 EDT 2012',
    'creator': 'Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G.',
    'section': 'Study',
    'title': 'Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific',
    'publisher': 'Systematics and Biodiversity',
    'identifier.study.tb1': None,
    'number': '',
    'identifier.study': '12713',
    'modificationDate': '2012-06-04',
    'historyNote': 'Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlDocumentWriter@645f9132 $Rev: 1060 $',
    'publicationDate': '2012',
    'contributor': 'Wiklund H.',
    'publicationName': 'Systematics and Biodiversity',
    'creationDate': '2012-05-09',
    'title.study': 'Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific',
    'subject': 'molecular phylogeny'}

Note that no attempt is made to prevent or account for key collision: :class:`~dendropy.datamodel.basemodel.Annotation` with the same name value will overwrite each other in the dictionary.
Custom control of the dictionary key/value generation can be specified via keyword arguments:

    ``key_attr``
        String specifying an Annotation object attribute name to be used
        as keys for the dictionary.

    ``key_func``
        Function that takes an Annotation object as an argument and returns
        the value to be used as a key for the dictionary.

    ``value_attr``
        String specifying an Annotation object attribute name to be used
        as values for the dictionary.

    ``value_func``
        Function that takes an Annotation object as an argument and returns
        the value to be used as a value for the dictionary.

For example::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    a = ds.annotations.values_as_dict(key_attr="prefixed_name")
    a = ds.annotations.values_as_dict(key_attr="prefixed_name", value_attr="namespace")
    a = ds.annotations.values_as_dict(key_func=lambda a: a.namespace + a.name)
    a = ds.annotations.values_as_dict(key_func=lambda a: a.namespace + a.name,
            value_attr="value")


As the collection returned by the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` method is an object of type :class:`~dendropy.datamodel.basemodel.AnnotationSet`, this can also be transformed to a dictionary.
For example::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    a = ds.annotations.findall(name_prefix="dc").values_as_dict()
    print a

will result in::

    {'publisher': 'Systematics and Biodiversity',
    'creator': 'Wiklund H., Altamira I.V., Glover A., Smith C., Baco A., & Dahlgren T.G.',
    'title': 'Systematics and biodiversity of Ophryotrocha (Annelida, Dorvilleidae) with descriptions of six new species from deep-sea whale-fall and wood-fall habitats in the north-east Pacific',
    'date': '2012-06-04',
    'contributor': 'Baco A.',
    'subject': 'molecular phylogeny'}

Note how only one entry for "contributor" is present: the others were overwritten/replaced.

Adding to, deleting, or modifying either the keys or the values of the dictionary returned by :meth:`~dendropy.datamodel.basemodel.Annotation.values_as_dict` in *no way* changes any of the original metadata: it is serves as snapshot copy of literal values of the metadata.


Deleting or Removing Metadata Annotations
-----------------------------------------

The :meth:`~dendropy.datamodel.basemodel.AnnotationSet.drop` method of :class:`~dendropy.datamodel.basemodel.AnnotationSet` objects takes search criteria similar to :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall`, but instead of returning the matched  :class:`~dendropy.datamodel.basemodel.Annotation` objects, it *removes* them from the parent collection.
For example, the following removes all metadata annotations with the name prefix "dc" from the |DataSet| object ``ds``::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    print "Original: %d items" % len(ds.annotations)
    removed = ds.annotations.drop(name_prefix="dc")
    print "Removed: %d items" % len(removed)
    print "Current: %d items" % len(ds.annotations)

and results in::

    Original: 30 items
    Removed: 16 items
    Current: 14 items

As can be seen, the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.drop` method returns the individual :class:`~dendropy.datamodel.basemodel.Annotation` removed as a new :class:`~dendropy.datamodel.basemodel.AnnotationSet` collection.
This is useful if you still want to use the removed :class:`~dendropy.datamodel.basemodel.Annotation` objects elsewhere.

As with the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` method, multiple keyword criteria can be specified::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
            "nexml")
    ds.annotations.drop(name_prefix="dc", name="contributor")

In addition, again similar in behavior to the :meth:`~dendropy.datamodel.basemodel.AnnotationSet.findall` method, *no* keyword arguments result in *all* the annotations being removed.
Thus, the following results in all metadata annotations being deleted from the |DataSet| object ``ds``::

    import dendropy
    ds = dendropy.DataSet.get_from_path("sample1.xml",
        "nexml")
    print "Original: %d items" % len(ds.annotations)
    removed = ds.annotations.drop()
    print "Removed: %d items" % len(removed)
    print "Current: %d items" % len(ds.annotations)

and results in::

    Original: 30 items
    Removed: 30 items
    Current: 0 items

Writing or Saving Metadata
==========================

When writing to NeXML format, all metadata annotations are preserved and can be fully round-tripped.
Currently, this is the only data format that allows for robust treatment of metadata.

Due to the fundamental limitations of the NEXUS/Newick format, metadata handling in this format is limited and rather idiosyncratic.
Currently, metadata will be written out as name-value pairs (separated by "=") in ampersand-prepended comments associated with the particular phylogenetic data object.
This syntax corresponds to the BEAST or FigTree style of metadata annotation.
However, this association might not be preserved.
For example, metadata annotations associated with edges and nodes of trees will be written out fully in NEXUS and NEWICK formats, but when read in again will all be associated with nodes.
The keyword argument ``annotations_as_nhx=True`` passed to the call to write the data in NEXUS/NEWICK format will result in a double ampersand prefix to the comment, thus (partially) conforming to NHX specifications.
Metadata associated with |DataSet| objects will be written in out in the same BEAST/FigTree/NHX syntax at the top of the file, while metadata associated with |TaxonNamespace| and |Taxon| objects will be written out immediately after the start of the Taxa Block and taxon labels respectively.
This is very fragile: for example, a metadata annotation *before* a taxon label will be associated with the *previous* taxon when being read in again.
As noted above, if metadata annotations are important for yourself, your workflow, or your task, then the NeXML format should be used rather than NEXUS or NEWICK.

