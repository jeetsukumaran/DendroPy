#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Infrastructure for object serialization and management.
"""

import re
import copy
from dendropy.utility import bibtex
from dendropy.utility import containers
from dendropy.utility import error

class DataObject(object):
    """
    Base class from which all classes that need to persist object attributes
    or other information as metadata.
    """

    def __init__(self, label=None, oid=None):
        self.label = label
        if oid is not None:
            self._oid = oid
        else:
            self._oid = self._default_oid()
        self._annotations = AnnotationSet(self)

    def _default_oid(self):
        "Returns default oid."
        return "d" + str(id(self))
    default_oid = property(_default_oid)

    def _get_oid(self):
        "Returns id."
        return self._oid
    def _set_oid(self, oid):
        """
        Sets oid to oid if oid is not None (normalized to conform
        to xs:NCName specs if neccessary), otherwise sets to some
        other oidue.
        """
        if oid is not None:
            self._oid = oid
        else:
            self._oid = self._default_oid()
    oid = property(_get_oid, _set_oid)

    def reset_oid(self):
        self._oid = self._default_oid
        return self._oid

    def __deepcopy__(self, memo):
        o = self.__class__(label=self.label, oid=None)
        memo[id(self)] = o
        memo[id(self._oid)] = o._oid
        memo[id(self.label)] = o.label
        for k, v in self.__dict__.iteritems():
            o.__dict__[copy.deepcopy(k, memo)] = copy.deepcopy(v, memo)
        return o

class AnnotatedDataObject(DataObject):
    """
    Base class from which all classes that need to persist object attributes
    or other information as metadata.
    """

    def __init__(self, label=None, oid=None):
        DataObject.__init__(self, label=label, oid=oid)
        self._annotations = AnnotationSet(self)

    def _get_annotations(self):
        if not hasattr(self, "_annotations"):
            self._annotations = AnnotationSet(self)
        return self._annotations
    def _set_annotations(self, annotations):
        if hasattr(self, "_annotations") \
                and annotations is self._annotations \
                and self._annotations.target is self:
            return
        if not isinstance(annotations, AnnotationSet):
            raise ValueError("Cannot set `annotations` to object of type `%s`" % type(annotations))
        old_target = annotations.target
        self._annotations = annotations
        self._annotations.target = self
        for a in self._annotations:
            if a.is_attribute and a._value[0] is old_target:
                a.target = self
    annotations = property(_get_annotations, _set_annotations)

    def __str__(self):
        return str(self.oid)

    def __deepcopy__(self, memo):

        # temporary disable deepcopying of annotations
        # until target object is created
        memo[id(self.annotations)] = None

        # copy object
        o = DataObject.__deepcopy__(self, memo)

        # reassign annotations
        # del memo[id(self.annotations)]
        # annotations_copy = copy.deepcopy(self.annotations, memo)
        # o._annotations = annotations_copy
        # memo[id(self.annotations)] = o._annotations
        del memo[id(self._annotations)]
        annotations_copy = copy.deepcopy(self._annotations, memo)
        o.annotations = annotations_copy
        memo[id(self._annotations)] = o._annotations

        # return
        return o

class Annotation(AnnotatedDataObject):
    """
    Metadata storage, composition and persistance, with the following attributes:

        - `name`
        - `value`
        - `datatype_hint`
        - `name_prefix`
        - `namespace`
        - `annotate_as_reference`
        - `is_hidden`

    """

    def parse_prefixed_name(prefixed_name, sep=":"):
        if sep not in prefixed_name:
            raise ValueError("'%s' is not a valid CURIE-standard qualified name" % prefixed_name)
        return prefixed_name.split(":", 1)
    parse_prefixed_name = staticmethod(parse_prefixed_name)

    def __init__(self,
            name,
            value,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_prefixed=False,
            is_attribute=False,
            annotate_as_reference=False,
            is_hidden=False,
            label=None,
            oid=None,
            ):
        AnnotatedDataObject.__init__(self, label=label, oid=oid)
        self._value = value
        self.is_attribute = is_attribute
        if name_is_prefixed:
            self.prefixed_name = name
            if name_prefix is not None:
                self._name_prefix = name_prefix
        else:
            self.name = name
            self._name_prefix = name_prefix
        self.datatype_hint = datatype_hint
        self._namespace = None
        self.namespace = namespace
        self.annotate_as_reference = annotate_as_reference
        self.is_hidden = is_hidden

    def __str__(self):
        return '%s="%s"' % (self.name, self.value)

    def __deepcopy__(self, memo):
        o = self.__class__(name=None, value=None)
        memo[id(self)] = o
        for k, v in self.__dict__.iteritems():
            if k not in ["label", "_oid"]:
                o.__dict__[copy.deepcopy(k, memo)] = copy.deepcopy(v, memo)
        # if o.is_attribute:
        #     o._value = (memo[i] for i in self._value)
        return o

    def is_match(self, **kwargs):
        match = True
        for k, v in kwargs.items():
            if k == "name_prefix":
                if self.name_prefix != v:
                    return False
            elif k == "prefixed_name":
                if self.prefixed_name != v:
                    return False
            elif k == "namespace":
                if self.namespace != v:
                    return False
            elif k == "value":
                if self.value != v:
                    return False
            elif hasattr(self, k):
                if getattr(self, k) != v:
                    return False
        return True

    def _get_value(self):
        if self.is_attribute:
            return getattr(*self._value)
        else:
            return self._value
    def _set_value(self, value):
        self._value = value
    value = property(_get_value, _set_value)

    def _get_name_prefix(self):
        if self._name_prefix is None:
            self._name_prefix = "dendropy"
        return self._name_prefix
    def _set_name_prefix(self, prefix):
        self._name_prefix = prefix
    name_prefix = property(_get_name_prefix, _set_name_prefix)

    def _get_namespace(self):
        if self._namespace is None:
            self._namespace = "http://packages.python.org/DendroPy/"
        return self._namespace
    def _set_namespace(self, prefix):
        self._namespace = prefix
    namespace = property(_get_namespace, _set_namespace)

    def _get_prefixed_name(self):
        return "%s:%s" % (self.name_prefix, self.name)
    def _set_prefixed_name(self, prefixed_name):
        self._name_prefix, self.name = Annotation.parse_prefixed_name(prefixed_name)
    prefixed_name = property(_get_prefixed_name, _set_prefixed_name)

class AnnotationSet(containers.OrderedSet):

    def __init__(self, target, *args):
        # for item in args:
        #     if not isinstance(item, Annotation):
        #         raise ValueError("Cannot add object of type `%s`" % type(item))
        containers.OrderedSet.__init__(self, *args)
        self.target = target

    def __str__(self):
        return "AnnotationSet([%s])" % ( ", ".join(str(a) for a in self))

    def __deepcopy__(self, memo):
        try:
            o = self.__class__(target=memo[id(self.target)])
        except KeyError:
            raise KeyError("deepcopy error: object id %d not found: %s" % (id(self.target), repr(self.target)))
        memo[id(self)] = o
        for a in self:
            o.add(copy.deepcopy(a, memo))
        return o

    def add_new(self,
            name,
            value,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_prefixed=False,
            is_attribute=False,
            annotate_as_reference=False,
            is_hidden=False):
        """
        Add an annotation, where:

            `name`
                The property/subject/field of the annotation (e.g. "color",
                "locality", "dc:citation")

            `value`
                The content of the annotation.

            `datatype_hint`
                Mainly for NeXML output (e.g. "xsd:string").

            `namespace_prefix`
                Mainly for NeXML output (e.g. "dc:").

            `namespace`
                Mainly for NeXML output (e.g. "http://www.w3.org/XML/1998/namespace").

            `name_is_prefixed`
                Mainly for NeXML *input*: name will be split into prefix and local part
                before storage (e.g., "dc:citations" will result in prefix = "dc" and
                name="citations")

            `is_attribute`
                If value is passed as a tuple of (object, "attribute_name") and this
                is True, then actual content will be the result of calling
                `getattr(object, "attribute_name")`.

            `annotate_as_reference`
                The value should be interpreted as a URI that points to content.

            `is_hidden`
                Do not write or print this annotation when writing data.

        """
        if not name_is_prefixed:
            if name_prefix is None and namespace is None:
                name_prefix = "dendropy"
                namespace = "http://packages.python.org/DendroPy/"
            elif name_prefix is None:
                raise TypeError("Cannot specify 'name_prefix' for unqualified name without specifying 'namespace'")
            elif namespace is None:
                raise TypeError("Cannot specify 'namespace' for unqualified name without specifying 'name_prefix'")
        else:
            if namespace is None:
                raise TypeError("Cannot specify qualified name without specifying 'namespace'")
        annote = Annotation(
                name=name,
                value=value,
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_prefixed=name_is_prefixed,
                is_attribute=is_attribute,
                annotate_as_reference=annotate_as_reference,
                is_hidden=is_hidden,
                )
        return self.add(annote)

    # def update(self, other):
    #     for item in other:
    #         if not isinstance(item, Annotation):
    #             raise ValueError("Cannot add object of type `%s`" % type(item))
    #     containers.OrderedSet.update(self, other)

    # def extend(self, other):
    #     for item in other:
    #         if not isinstance(item, Annotation):
    #             raise ValueError("Cannot add object of type `%s`" % type(item))
    #     containers.OrderedSet.extend(self, other)

    # def add(self, item):
    #     if not isinstance(item, Annotation):
    #         raise ValueError("Cannot add object of type `%s`" % type(item))
    #     containers.OrderedSet.add(self, item)

    # def append(self, item):
    #     if not isinstance(item, Annotation):
    #         raise ValueError("Cannot add object of type `%s`" % type(item))
    #     containers.OrderedSet.append(self, item)

    def add_bound_attribute(self,
            attr_name,
            annotation_name=None,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_prefixed=False,
            annotate_as_reference=False,
            is_hidden=False,
            owner_instance=None,
            ):
        """
        Add an attribute of an object as an annotation. The value of the
        annotation will be dynamically bound to the value of the attribute.

            `attr_name`
                The (string) name of the attribute to be used as the source of the
                content or value of the annotation.

            `annotation_name`
                Use this string as the annotation field/name rather than the attribute
                name.

            `datatype_hint`
                Mainly for NeXML output (e.g. "xsd:string").

            `namespace_prefix`
                Mainly for NeXML output (e.g. "dc:").

            `namespace`
                Mainly for NeXML output (e.g. "http://www.w3.org/XML/1998/namespace").

            `name_is_prefixed`
                Mainly for NeXML *input*: name will be split into prefix and local part
                before storage (e.g., "dc:citations" will result in prefix = "dc" and
                name="citations")

            `annotate_as_reference`
                The value should be interpreted as a URI that points to content.

            `is_hidden`
                Do not write or print this annotation when writing data.

            `owner_instance`
                The object whose attribute is to be used as the value of the
                annotation. Defaults to `self.target`.

        """
        if annotation_name is None:
            annotation_name = attr_name
        if owner_instance is None:
            owner_instance = self.target
        if not hasattr(owner_instance, attr_name):
            raise AttributeError(attr_name)
        if not name_is_prefixed:
            if name_prefix is None and namespace is None:
                name_prefix = "dendropy"
                namespace = "http://packages.python.org/DendroPy/"
            elif name_prefix is None:
                raise TypeError("Cannot specify 'name_prefix' for unqualified name without specifying 'namespace'")
            elif namespace is None:
                raise TypeError("Cannot specify 'namespace' for unqualified name without specifying 'name_prefix'")
        else:
            if namespace is None:
                raise TypeError("Cannot specify qualified name without specifying 'namespace'")
        annote = Annotation(
                name=annotation_name,
                value=(owner_instance, attr_name),
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_prefixed=name_is_prefixed,
                is_attribute=True,
                annotate_as_reference=annotate_as_reference,
                is_hidden=is_hidden,
                )
        return self.add(annote)

    def add_citation(self,
            citation,
            read_as="bibtex",
            store_as="bibtex",
            name_prefix=None,
            namespace=None,
            is_hidden=False):
        """
        Add a citation as an annotation.

            `citation`
                This can be a string representing a BibTex entry, a dictionary
                with BibTex fields as keys and contents as values, or a
                dendropy.utility.BibTex.BibTexEntry object.

            `read_as`
                Specifies the format/schema/structure of the citation. Currently
                only supports 'bibtex'.

            `store_as`
                Specifies how to record the citation, takes one of the
                following strings as values:

                    "bibtex"
                        A set of annotations, where each BibTex field becomes a
                        separate annotation.

                    "prism"
                        A set of of PRISM (Publishing Requirements for Industry
                        Standard Metadata) annotations.

                    "dublin"
                        A set of of Dublic Core annotations.

                Defaults to "bibtex".

            `name_prefix`
                Mainly for NeXML output (e.g. "dc:").

            `namespace`
                Mainly for NeXML output (e.g. "http://www.w3.org/XML/1998/namespace").

            `is_hidden`
                Do not write or print this annotation when writing data.

        """
        if read_as == "bibtex":
            return self.add_bibtex(citation=citation,
                    store_as=store_as,
                    name_prefix=name_prefix,
                    namespace=namespace,
                    is_hidden=is_hidden)
        else:
            raise ValueError("Source format '%s' is not supported" % read_as)

    def add_bibtex(self,
            citation,
            store_as="bibtex",
            name_prefix=None,
            namespace=None,
            is_hidden=False):
        """
        Add a citation as an annotation.

            `citation`
                This can be a string representing a BibTex entry, a dictionary
                with BibTex fields as keys and contents as values, or a
                dendropy.utility.BibTex.BibTexEntry object.

            `store_as`
                Specifies how to record the citation:

                    "bibtex"
                        A set of annotations, where each BibTex field becomes a
                        separate annotation.

                    "prism"
                        A set of of PRISM (Publishing Requirements for Industry
                        Standard Metadata) annotations.

                    "dublin"
                        A set of of Dublic Core annotations.

                Defaults to "bibtex".

            `name_prefix`
                Mainly for NeXML output (e.g. "dc:").

            `namespace`
                Mainly for NeXML output (e.g. "http://www.w3.org/XML/1998/namespace").

            `is_hidden`
                Do not write or print this annotation when writing data.

        """
        bt = bibtex.BibTexEntry(citation)
        bt_dict = bt.fields_as_dict()

        if name_prefix is None and namespace is not None:
            raise TypeError("Cannot specify 'name_prefix' for unqualified name without specifying 'namespace'")
        elif namespace is None and name_prefix is not None:
            raise TypeError("Cannot specify 'namespace' for unqualified name without specifying 'name_prefix'")

        if store_as.lower().startswith("bibtex"):
            if name_prefix is None and namespace is None:
                name_prefix = "bibtex"
                namespace = "http://www.edutella.org/bibtex#"
            self.add_new(
                    name="bibtype",
                    value=bt.bibtype,
                    datatype_hint="xsd:string",
                    name_prefix=name_prefix,
                    namespace=namespace,
                    name_is_prefixed=False,
                    is_attribute=False,
                    annotate_as_reference=False,
                    is_hidden=is_hidden)
            self.add_new(
                    name="citekey",
                    value=bt.citekey,
                    datatype_hint="xsd:string",
                    name_prefix=name_prefix,
                    namespace=namespace,
                    name_is_prefixed=False,
                    is_attribute=False,
                    annotate_as_reference=False,
                    is_hidden=is_hidden)
            for entry_key, entry_value in bt_dict.items():
                self.add_new(
                        name=entry_key,
                        value=entry_value,
                        datatype_hint="xsd:string",
                        name_prefix=name_prefix,
                        namespace=namespace,
                        name_is_prefixed=False,
                        is_attribute=False,
                        annotate_as_reference=False,
                        is_hidden=is_hidden)
        # elif store_as.lower().startswith("bibtex-record"):
        #     if name_prefix is None and namespace is None:
        #         name_prefix = "dendropy"
        #         namespace = "http://packages.python.org/DendroPy/"
        #     self.add_new(
        #             name="bibtex",
        #             value=bt.as_compact_bibtex(),
        #             datatype_hint="xsd:string",
        #             name_is_prefixed=False,
        #             name_prefix=name_prefix,
        #             namespace=namespace,
        #             is_attribute=False,
        #             annotate_as_reference=False,
        #             is_hidden=is_hidden)
        elif store_as.lower().startswith("prism"):
            prism_map = {
                    'volume': bt_dict.get('volume', None),
                    'publicationName':  bt_dict.get('journal', None),
                    'pageRange': bt_dict.get('pages', None),
                    'publicationDate': bt_dict.get('year', None),
                    }
            if name_prefix is None and namespace is None:
                name_prefix = "prism"
                namespace = "http://prismstandard.org/namespaces/1.2/basic/"
            for field, value in prism_map.items():
                if value is None:
                    continue
                self.add_new(
                        name=field,
                        value=value,
                        datatype_hint="xsd:string",
                        name_prefix=name_prefix,
                        namespace=namespace,
                        name_is_prefixed=False,
                        is_attribute=False,
                        annotate_as_reference=False,
                        is_hidden=is_hidden)
        elif store_as.lower().startswith("dublin"):
            dc_map = {
                    'title': bt_dict.get('title', None),
                    'creator':  bt_dict.get('author', None),
                    'publisher': bt_dict.get('journal', None),
                    'date': bt_dict.get('year', None),
                    }
            if name_prefix is None and namespace is None:
                name_prefix = "dc"
                namespace = "http://purl.org/dc/elements/1.1/"
            for field, value in dc_map.items():
                if value is None:
                    continue
                self.add_new(
                        name=field,
                        value=value,
                        datatype_hint="xsd:string",
                        name_is_prefixed=False,
                        name_prefix=name_prefix,
                        namespace=namespace,
                        is_attribute=False,
                        annotate_as_reference=False,
                        is_hidden=is_hidden)
        else:
            raise ValueError("Unrecognized composition specification: '%s'" % store_as)

    def findall(self, **kwargs):
        """
        Returns AnnotationSet of Annotation objects associated with self.target
        that match based on *all* criteria specified in keyword arguments::

            >>> notes = tree.annotations.findall(name="color")
            >>> notes = tree.annotations.findall(namespace="http://packages.python.org/DendroPy/")
            >>> notes = tree.annotations.findall(namespace="http://packages.python.org/DendroPy/",
                                          name="color")
            >>> notes = tree.annotations.findall(name_prefix="dc")
            >>> notes = tree.annotations.findall(prefixed_name="dc:color")

        If no matches are found, the return AnnotationSet is empty.

        If no keyword arguments are given, *all* annotations are returned::

            >>> notes = tree.annotations.findall()

        """
        results = []
        for a in self:
            if a.is_match(**kwargs):
                results.append(a)
        results = AnnotationSet(self.target, results)
        return results

    def find(self, **kwargs):
        """
        Returns the *first* Annotation associated with self.target
        which matches based on *all* criteria specified in keyword arguments::

            >>> note = tree.annotations.find(name="color")
            >>> note = tree.annotations.find(name_prefix="dc", name="color")
            >>> note = tree.annotations.find(prefixed_name="dc:color")

        If no match is found, None is returned.

        If no keyword arguments are given, a TypeError is raised.
        """
        if "default" in kwargs:
            default = kwargs["default"]
            del kwargs["default"]
        else:
            default = None
        if not kwargs:
            raise TypeError("Search criteria not specified")
        for a in self:
            if a.is_match(**kwargs):
                return a
        return default

    def get_value(self, name, default=None):
        """
        Returns the *value* of the *first* Annotation associated with
        self.target which has `name` in the name field.

        If no match is found, then `default` is returned.
        """
        for a in self:
            if a.is_match(name=name):
                return a.value
        return default

    def drop(self, **kwargs):
        """
        Removes Annotation objects that match based on *all* criteria specified
        in keyword arguments.

        Remove all annotation objects with `name` ==
        "color"::

            >>> tree.annotations.drop(name="color")

        Remove all annotation objects with `namespace` ==
        "http://packages.python.org/DendroPy/"::

            >>> tree.annotations.drop(namespace="http://packages.python.org/DendroPy/")

        Remove all annotation objects with `namespace` ==
        "http://packages.python.org/DendroPy/" *and* `name` == "color"::

            >>> tree.annotations.drop(namespace="http://packages.python.org/DendroPy/",
                    name="color")

        Remove all annotation objects with `name_prefix` == "dc"::

            >>> tree.annotations.drop(name_prefix="dc")

        Remove all annotation objects with `prefixed_name` == "dc:color"::

            >>> tree.annotations.drop(prefixed_name="dc:color")

        If no keyword argument filter criteria are given, *all* annotations are
        removed::

            >>> tree.annotations.drop()

        """
        to_remove = []
        for a in self:
            if a.is_match(**kwargs):
                to_remove.append(a)
        for a in to_remove:
            self.remove(a)
        return AnnotationSet(self.target, to_remove)

    def values_as_dict(self, **kwargs):
        """
        Returns annotation set as a dictionary. The keys and values for the dictionary will
        be generated based on the following keyword arguments:

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

        At most one of ``key_attr`` or ``key_func`` can be specified. If neither
        is specified, then by default the keys are generated from Annotation.name.
        At most one of ``value_attr`` or ``value_func`` can be specified. If neither
        is specified, then by default the values are generated from Annotation.value.
        Key collisions will result in the dictionary entry for that key being
        overwritten.
        """
        if "key_attr" in kwargs and "key_func" in kwargs:
            raise TypeError("Cannot specify both 'key_attr' and 'key_func'")
        elif "key_attr" in kwargs:
            key_attr = kwargs["key_attr"]
            key_func = lambda a: getattr(a, key_attr)
        elif "key_func" in kwargs:
            key_func = kwargs["key_func"]
        else:
            key_func = lambda a: a.name
        if "value_attr" in kwargs and "value_func" in kwargs:
            raise TypeError("Cannot specify both 'value_attr' and 'value_func'")
        elif "value_attr" in kwargs:
            value_attr = kwargs["value_attr"]
            value_func = lambda a: getattr(a, value_attr)
        elif "value_func" in kwargs:
            value_func = kwargs["value_func"]
        else:
            value_func = lambda a: a.value
        d = {}
        for a in self:
            d[key_func(a)] = value_func(a)
        return d
