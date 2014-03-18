#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Infrastructure for phylogenetic data objects.
"""

import os
import copy
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.utility import container
from dendropy.utility import bibtex

##############################################################################
## DataObject

class DataObject(object):

    """
    Base class for all phylogenetic data objects.
    """

    def __init__(self, label=None):
        self._label=label

    def _get_label(self):
        return self._label
    def _set_label(self, v):
        self._label = v
    label = property(_get_label, _set_label)

   # def __copy__(self):
   #      cls = self.__class__
   #      result = cls.__new__(cls)
   #      result.__dict__.update(self.__dict__)
   #      return result

    # def __deepcopy__(self, memo=None):
    #     # cls = self.__class__
    #     # result = cls.__new__(cls)
    #     # memo[id(self)] = result
    #     # for k, v in self.__dict__.items():
    #     #     setattr(result, k, deepcopy(v, memo))
    #     # return result
    #     if memo is None:
    #         memo = {}
    #     try:
    #         o = memo[id(self)]
    #     except KeyError:
    #         # o = type(self).__new__(self.__class__)
    #         o = self.__class__.__new__(self.__class__)
    #         memo[id(self)] = o
    #     for k in self.__dict__:
    #         # o.__dict__[copy.deepcopy(k, memo)] = copy.deepcopy(self.__dict__[k], memo)
    #         o.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
    #         # setattr(o, k, copy.deepcopy(self.__dict__[k], memo))
    #         # setattr(o, k, copy.deepcopy(self.__dict__[k], memo))
    #     return o

##############################################################################
## Readable

class Readable(object):
    """
    Mixin class which all classes that require deserialization should subclass.
    """

    def parse_from_stream(cls, stream, schema, **kwargs):
        """
        Subclasses need to implement this method to create
        and return and instance of themselves read from the
        stream.
        """
        raise NotImplementedError
    parse_from_stream = classmethod(parse_from_stream)

    def get_from_stream(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from file-like
        object `src`.
        """
        return cls.parse_from_stream(stream=src,
                schema=schema,
                **kwargs)
    get_from_stream = classmethod(get_from_stream)

    def get_from_path(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from file
        specified by string `src`.
        """
        fsrc = open(src, "rU")
        return cls.parse_from_stream(stream=fsrc,
                schema=schema,
                **kwargs)
    get_from_path = classmethod(get_from_path)

    def get_from_string(cls, src, schema, **kwargs):
        """
        Factory method to return new object of this class from string `src`.
        """
        ssrc = StringIO(src)
        return cls.parse_from_stream(stream=ssrc,
                schema=schema,
                **kwargs)
    get_from_string = classmethod(get_from_string)

    def get_from_url(cls, src, schema, strip_markup=False, **kwargs):
        """
        Factory method to return a new object of this class from
        URL given by `src`.
        """
        text = read_url(src, strip_markup=strip_markup)
        ssrc = StringIO(text)
        try:
            return cls.parse_from_stream(stream=ssrc,
                    schema=schema,
                    **kwargs)
        except error.DataParseError:
            sys.stderr.write(text)
            raise
    get_from_url = classmethod(get_from_url)

    def __init__(self, *args, **kwargs):
        pass

    def process_source_kwargs(self, **kwargs):
        """
        If `stream` is specified, then process_source_kwargs:

            # checks that `schema` keyword argument is specified, and then
            # calls self.read()

        If `stream` is not specified, then nothing happens: unless the
        data object was populated through other means, and empty data
        object will the result (typically used as a starting point
        for population unit by unit).
        """
        if "stream" in kwargs:
            stream = kwargs["stream"]
            del(kwargs["stream"])
            schema = require_format_from_kwargs(kwargs)
            self.read(stream=stream, schema=schema, **kwargs)
        # else:
        #     from pudb import set_trace; set_trace()
        # elif "source_string" in kwargs:
        #     as_str = kwargs["source_string"]
        #     del(kwargs["source_string"])
        #     kwargs["stream"] = StringIO(as_str)
        #     return self.process_source_kwargs(**kwargs)
        # elif "source_file" in kwargs:
        #     fp = kwargs["source_filepath"]
        #     fo = open(os.path.expandvars(os.path.expanduser(filepath)), "rU")
        #     del(kwargs["source_filepath"])
        #     kwargs["stream"] = fo
        #     return self.process_source_kwargs(**kwargs)

    def read(self, stream, schema, **kwargs):
        """
        Populates/constructs objects of this type from `schema`-formatted
        data in the file-like object source `stream`.
        """
        raise NotImplementedError

    def read_from_stream(self, fileobj, schema, **kwargs):
        """
        Reads from file (exactly equivalent to just `read()`, provided
        here as a separate method for completeness.
        """
        return self.read(stream=fileobj, schema=schema, **kwargs)

    def read_from_path(self, filepath, schema, **kwargs):
        """
        Reads from file specified by `filepath`.
        """
        f = open(os.path.expandvars(os.path.expanduser(filepath)), "rU")
        return self.read(stream=f, schema=schema, **kwargs)

    def read_from_string(self, src_str, schema, **kwargs):
        """
        Reads a string object.
        """
        s = StringIO(src_str)
        return self.read(stream=s, schema=schema, **kwargs)

    def read_from_url(self, url, schema, **kwargs):
        """
        Reads a URL source.
        """
        src_str = read_url(url)
        s = StringIO(src_str)
        return self.read(stream=s, schema=schema, **kwargs)

##############################################################################
## Writeable

class Writeable(object):
    """
    Mixin class which all classes that require serialization should subclass.
    """

    def __init__(self, *args, **kwargs):
        pass

    def write(self, stream, schema, **kwargs):
        """
        Writes the object to the file-like object `stream` in `schema`
        schema.
        """
        raise NotImplementedError

    def write_to_stream(self, dest, schema, **kwargs):
        """
        Writes to file-like object `dest`.
        """
        return self.write(stream=dest, schema=schema, **kwargs)

    def write_to_path(self, dest, schema, **kwargs):
        """
        Writes to file specified by `dest`.
        """
        f = open(os.path.expandvars(os.path.expanduser(dest)), "w")
        return self.write(stream=f, schema=schema, **kwargs)

    def as_string(self, schema, **kwargs):
        """
        Composes and returns string representation of the data.
        """
        s = StringIO()
        self.write(stream=s, schema=schema, **kwargs)
        return s.getvalue()

##############################################################################
## Annotable

class Annotable(object):
    """
    Mixin class which all classes that need to persist object attributes
    or other information as metadata should subclass.
    """

    def __init__(self):
        pass

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
            raise ValueError("Cannot set 'annotations' to object of type '{}'".format(type(annotations)))
        old_target = annotations.target
        self._annotations = annotations
        self._annotations.target = self
        for a in self._annotations:
            if a.is_attribute and a._value[0] is old_target:
                a.target = self
    annotations = property(_get_annotations, _set_annotations)

    def __str__(self):
        return str(self.oid)

    def copy_annotations_from(self, other):
        if hasattr(other, "_annotations"):
            for annote in other._annotations:
                a2 = copy.deepcopy(annote)
                if a2.is_attribute and a2.value[0] is other:
                    a2.value[1] = self
                self.annotations.add(a2)

    def __deepcopy__(self, memo=None):
        # ensure clone map
        if memo is None:
            memo = {}
        # get or create clone of self
        try:
            o = memo[id(self)]
        except KeyError:
            # create object without initialization
            # o = type(self).__new__(self.__class__)
            o = self.__class__.__new__(self.__class__)
            # store
            memo[id(self)] = o
        # create annotations
        if hasattr(self, "_annotations") and id(self._annotations) not in memo:
            o._annotations = copy.deepcopy(self._annotations, memo)
            memo[id(self._annotations)] = o._annotations

        # return
        return o

##############################################################################
## Annotation

class Annotation(Annotable):
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
            raise ValueError("'{}' is not a valid CURIE-standard qualified name".format(prefixed_name))
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
            ):
        Annotable.__init__(self)
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
        return "{}='{}'".format(self.name, self.value)

    def __deepcopy__(self, memo):
        o = Annotable.__deepcopy__(self, memo)
        memo[id(self)] = o
        for k in self.__dict__:
            if k not in o.__dict__: # do not add attributes already added by base class
                o.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
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
        return "{}:{}".format((self.name_prefix, self.name))
    def _set_prefixed_name(self, prefixed_name):
        self._name_prefix, self.name = Annotation.parse_prefixed_name(prefixed_name)
    prefixed_name = property(_get_prefixed_name, _set_prefixed_name)

##############################################################################
## AnnotationSet

class AnnotationSet(container.OrderedSet):

    def __init__(self, target, *args):
        container.OrderedSet.__init__(self, *args)
        self.target = target

    def __str__(self):
        return "AnnotationSet([{}])".format(( ", ".join(str(a) for a in self)))

    def __deepcopy__(self, memo):
        try:
            o = self.__class__(target=memo[id(self.target)])
        except KeyError:
            raise KeyError("deepcopy error: object id {} not found: {}".format((id(self.target), repr(self.target))))
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
            raise ValueError("Source format '{}' is not supported".format(read_as))

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
            raise ValueError("Unrecognized composition specification: '{}'".format(store_as))

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

