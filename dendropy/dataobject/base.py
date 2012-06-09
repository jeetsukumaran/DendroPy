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

class DataObject(object):
    pass

class Annotated(DataObject):
    """
    Base class from which all classes that need to persist object attributes
    beyond the core elements (such as id, label, etc.) will derive.
    """

    def __init__(self, label=None, oid=None):
        self.label = label
        if oid is not None:
            self._oid = oid
        else:
            self._oid = self._default_oid()
        self.annotations = set()

    def _default_oid(self):
        "Returns default oid."
        return "x" + str(id(self))

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

    def add_annotation(self,
            name,
            value,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_qualified=False,
            value_is_attribute=False):
        """
        Add an annotation to self, where:

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

            `name_is_qualified`
                Mainly for NeXML *input*: name will be split into prefix and local part
                before storage (e.g., "dc:citations" will result in prefix = "dc" and
                name="citations")

            `value_is_attribute`
                If value is passed as a tuple of (object, "attribute_name") and this
                is True, then actual content will be the result of calling
                `getattr(object, "attribute_name")`.

        """
        if not name_is_qualified:
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
                name_is_qualified=name_is_qualified,
                value_is_attribute=value_is_attribute,
                )
        self.annotations.add(annote)
        return annote

    def add_attribute_annotation(self,
            attr_name,
            annotate_as=None,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_qualified=False,
            ):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation. The value of the annotation will
        be dynamically bound to the value of the attribute.

            `attr_name`
                The name of the attribute to be used as the source of the
                content or value of the annotation.

            `annotate_as`
                Use this string as the annotation field/name rather than the attribute
                name.

            `datatype_hint`
                Mainly for NeXML output (e.g. "xsd:string").

            `namespace_prefix`
                Mainly for NeXML output (e.g. "dc:").

            `namespace`
                Mainly for NeXML output (e.g. "http://www.w3.org/XML/1998/namespace").

            `name_is_qualified`
                Mainly for NeXML *input*: name will be split into prefix and local part
                before storage (e.g., "dc:citations" will result in prefix = "dc" and
                name="citations")

        """
        if annotate_as is None:
            annotate_as = attr_name
        if not hasattr(self, attr_name):
            raise AttributeError(attr_name)
        if not name_is_qualified:
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
                name=annotate_as,
                value=(self, attr_name),
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_qualified=name_is_qualified,
                value_is_attribute=True,
                )
        self.annotations.add(annote)
        return annote

    def has_annotations(self):
        """
        Returns True if there are attributes to be persisted as
        annotations.
        """
        return bool(len(self.annotations) > 0)

    def get_annotations(self, **kwargs):
        """
        Returns list of Annotation objects associated with self that match
        based on *all* criteria specified in keyword arguments::

            >>> notes = a.get_annotations(name="color")
            >>> notes = a.get_annotations(namespace="http://packages.python.org/DendroPy/")
            >>> notes = a.get_annotations(namespace="http://packages.python.org/DendroPy/",
                                          name="color")
            >>> notes = a.get_annotations(name_prefix="dc")
            >>> notes = a.get_annotations(qualified_name="dc:color")

        If no keyword arguments are given, *all* annotations are returned::

            >>> notes = a.get_annotations()

        """
        results = []
        for a in self.annotations:
            if a.is_match(**kwargs):
                results.append(a)
        return results

    def remove_annotations(self, **kwargs):
        """
        Removes Annotation objects associated with self that match
        based on *all* criteria specified in keyword arguments.

        Remove all annotation objects associated with self with `name` ==
        "color"::

            >>> a.remove_annotations(name="color")

        Remove all annotation objects associated with self with `namespace` ==
        "http://packages.python.org/DendroPy/"::

            >>> a.remove_annotations(namespace="http://packages.python.org/DendroPy/")

        Remove all annotation objects associated with self with `namespace` ==
        "http://packages.python.org/DendroPy/" *and* `name` == "color"::

            >>> a.remove_annotations(namespace="http://packages.python.org/DendroPy/",
                    name="color")

        Remove all annotation objects associated with self with `name_prefix` ==
        "dc"::

            >>> a.remove_annotations(name_prefix="dc")

        Remove all annotation objects associated with self with `qualified_name` ==
        "dc:color"::

            >>> a.remove_annotations(qualified_name="dc:color")

        If no keyword argument filter criteria are given, *all* annotations are
        removed::

            >>> a.remove_annotations()

        """
        to_remove = []
        for a in self.annotations:
            if a.is_match(**kwargs):
                to_remove.append(a)
        for a in to_remove:
            self.annotations.remove(a)

class Annotation(DataObject):
    "Tracks the basic information need to serialize an attribute correctly."

    def parse_qualified_name(qualified_name, sep=":"):
        if sep not in qualified_name:
            raise ValueError("'%s' is not a valid CURIE-standard qualified name" % qualified_name)
        return qualified_name.split(":", 1)
    parse_qualified_name = staticmethod(parse_qualified_name)

    def __init__(self,
            name,
            value,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_qualified=False,
            value_is_attribute=False,
            ):
        DataObject.__init__(self)
        self._value = value
        self.value_is_attribute = value_is_attribute
        if name_is_qualified:
            self.qualified_name = name
            if name_prefix is not None:
                self._name_prefix = name_prefix
        else:
            self.name = name
            self._name_prefix = name_prefix
        self.datatype_hint = datatype_hint
        self._namespace = None
        self.namespace = namespace

    def is_match(self, **kwargs):
        match = True
        for k, v in kwargs.items():
            if k == "name_prefix":
                if self.name_prefix != v:
                    return False
            elif k == "qualified_name":
                if self.qualified_name != v:
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
        if self.value_is_attribute:
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

    def _get_qualified_name(self):
        return "%s:%s" % (self.name_prefix, self.name)
    def _set_qualified_name(self, qualified_name):
        self._name_prefix, self.name = Annotation.parse_qualified_name(qualified_name)
    qualified_name = property(_get_qualified_name, _set_qualified_name)

