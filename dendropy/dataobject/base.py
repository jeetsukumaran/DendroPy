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
    """
    Base class for elements that will be serialized.
    """
    def __init__(self):
        self.attributes = []    ## extra attributes to be written during XML serialization
        # self.extensions = []    ## DOM extensions to be written during XML serialization

class Annotated(DataObject):
    """
    Base class from which all classes that need to persist object attributes
    beyond the core elements (such as id, label, etc.) will derive.
    """

    def __init__(self):
        pass

    def _create(self):
        DataObject.__init__(self)
        self.annotations = set()

    def __getattr__(self, name):
        if name == "annotations":
            self._create()
            return self.annotations
        raise AttributeError(name)

    def add_annotation(self,
            name,
            value,
            datatype_hint=None,
            name_prefix=None,
            namespace=None,
            name_is_qualified=False):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        annote = Annotation(
                name=name,
                value=value,
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_qualified=name_is_qualified,
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
        persisted as an annotation.
        """
        if annotate_as is None:
            annotate_as = attr_name
        if not hasattr(self, attr_name):
            raise AttributeError(attr_name)
        value = getattr(self, attr_name)
        return self.add_annotation(
                name=annotate_as,
                value=value,
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_qualified=name_is_qualified,
                )

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

            >>> a.get_annotations(name="color")
            >>> a.get_annotations(namespace="http://packages.python.org/DendroPy/")
            >>> a.get_annotations(namespace="http://packages.python.org/DendroPy/",
                    name="color")
            >>> a.get_annotations(name_prefix="dc")
            >>> a.get_annotations(qualified_name="dc:color")

        """
        results = []
        for a in self.annotations:
            if a.is_match(**kwargs):
                results.append(a)
        return results

    def remove_annotations(self, **kwargs):
        """
        Removes Annotation objects associated with self that match
        based on *all* criteria specified in keyword arguments::

            >>> a.remove_annotations(name="color")
            >>> a.remove_annotations(namespace="http://packages.python.org/DendroPy/")
            >>> a.remove_annotations(namespace="http://packages.python.org/DendroPy/",
                    name="color")
            >>> a.remove_annotations(name_prefix="dc")
            >>> a.remove_annotations(qualified_name="dc:color")

        """
        for a in self.annotations:
            if a.is_match(**kwargs):
                self.annotations.remove(a)

class Annotation(Annotated):
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
            name_is_qualified=False
            ):
        self.value = value
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
                if self.name_prefix == v:
                    return False
            elif k == "qualified_name":
                if self.qualified_name == v:
                    return False
            elif k == "namespace":
                if self.namespace == v:
                    return False
            elif hasattr(self, k):
                if getattr(self, k) != v:
                    return False
        return True

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

class Labelled(Annotated):
    "Provides for getting and setting of an object label."

    def __init__(self, label=None):
        """
        __init__ calls Annotated.__init__, and then, if keyword
        argument `label` is given, assigns it to self.label.
        """
        self.label = label

class IdTagged(Labelled):
    """
    Provides infrastructure for the maintenance of a unique object id,
    ensuring that this will never be None.
    """

    instances = 0

    # def normalize_id(id_str):
    #     """
    #     Given a string `id_str`, this returns a xs:NCName compliant
    #     version of the string: (Letter | '_' | ':')
    #     (NameChar)*. NameChar is given by : Letter | Digit | '.' | '-'
    #     | '_' | ':'
    #     """
    #     if len(id_str) > 0:
    #         f = id_str[0]
    #     else:
    #         f = '_' + str(id(id_str))
    #     if not (f.isalpha or f == '_' or f == ':'):
    #         id_str = '_' + id_str
    #     id_str = re.sub('[^\w\d\-\.]', '', id_str)
    #     return id_str

    # normalize_id = staticmethod(normalize_id)

    def __init__(self, label=None, oid=None, **kwargs):
        """
        __init__ calls Labelled.__init__, and assigns element id if
        given.
        """
        self.label = label
        IdTagged.instances += 1
        if oid is not None:
            self._oid = oid
        else:
            self._oid = self._default_oid()

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

    def __str__(self):
        "String representation of the object: it's id."
        return str(self.oid)
