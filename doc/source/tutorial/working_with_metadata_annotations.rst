*************************************
Working with Metadata and Annotations
*************************************

|DendroPy| provides a rich infrastructure for decorating most types of phylogenetic objects (e.g., the |DataSet|, |TaxonSet|, |Taxon| |TreeList|, |Tree|, and various |CharacterMatrix| classes) with metadata information.
These phylogenetic objects have an attribute, :attr:`annotations`, that is an instance of the :class:`~dendropy.dataobject.base.AnnotationSet` class, which is an iterable (derived from :class:`set`) that serves to manage a collection of :class:`~dendropy.dataobject.base.Annotation` objects.
Each :class:`~dendropy.dataobject.base.Annotation` object tracks a single annotation element.
These annotations will be rendered as ``meta`` elements when writing to NeXML format or "hot comments" when writing to NEXUS/NEWICK format.
Note that full and robust expression of metadata annotations, including stable and consistent round-tripping of information, can only be achieved while in the NeXML format.

Metadata Annotation Creation
=============================

Direct Composition with Literal Values
--------------------------------------

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
For convenience, you can specify the name of the annotation with the name prefix prepended by specifying "``name_is_qualified=True``", thought the namespace must still be provided separately::

    >>> tree.annotations.add_new(
    ... name="dc:subject",
    ... value="Python phylogenetics",
    ... name_is_qualified=True,
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

Sometimes, you may want to annotate an object with metadata, but do not want it to be printed or written out.
Passing the ``is_hidden=True`` argument will result in the annotation being suppressed in all output::

    >>> tree.annotations.add_new(
    ... name="subject",
    ... value="Python phylogenetics",
    ... name_prefix="dc",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... is_hidden=True,
    ... )

The ``is_hidden`` attribute of the an :class:`~dendropy.dataobject.base.Annotation` object can also be set directly::

    >>> subject_annotations = tree.annotations.get(name="citation")
    >>> for a in subject_annotations:
    ...    a.is_hidden = True

Dynamically Binding Annotation Values to Object Attribute Values
----------------------------------------------------------------

In some cases, instead of "hard-wiring" in metadata for an object, you may want to write out metadata that takes its value from the value of an attribute of the object.
The :meth:`~dendropy.dataobject.base.Annotation.add_bound_attribute` method allows you to do this.
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


By default, the :meth:`~dendropy.dataobject.base.Annotation.add_bound_attribute` method uses the name of the attribute as the name of the annotation.
The "``annotate_as``" argument allows you explictly set the name of the annotation.
In addition, the method call also supports the other customization arguments of the :meth:`~dendropy.dataobject.base.Annotation.add_new` method: "``datatype_hint``", "``name_prefix``", "``namespace``", "``name_is_qualified``", "``compose_as_reference``", "``is_hidden``", etc.::

    >>> tree.source_uri = None
    >>> tree.annotations.add_bound_attribute(
    ... "source_uri",
    ... annotate_as="dc:subject",
    ... namespace="http://purl.org/dc/elements/1.1/",
    ... compose_as_reference=True)

Adding Citation Metadata
------------------------

You can add citation annotations using the :meth:`~dendropy.dataobject.base.Annotation.add_citation` method.
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
        This is the default.

    - "``prism``"
        A subset of the BibTex fields gets recorded as a set of PRISM (Publishing Requirements for Industry Standard Metadata) annotations, one per field.


    - "``dublin``"
        A subset of the BibTex fields gets recorded as a set of Dublin Core (Publishing Requirements for Industry Standard Metadata) annotations, one per field.

In addition, the method call also supports some of the other customization arguments of the :meth:`~dendropy.dataobject.base.Annotation.add_new` method:  "``name_prefix``", "``namespace``", "``name_is_qualified``", "``is_hidden``".

.. literalinclude:: /examples/bibtex_annotations3.py
