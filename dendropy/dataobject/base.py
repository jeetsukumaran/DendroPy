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
        self.annotations = []

    def store_annotation(self,
            key,
            value,
            datatype_hint=None,
            namespace_map=None,
            namespace_key=None):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        annote = Annotation(
                key=key,
                value=value,
                datatype_hint=datatype_hint,
                namespace_map=namespace_map,
                namespace_key=namespace_key,
                )
        self.annotations.append(annote)
        return annote

    def __getattr__(self, name):
        if name == "annotations":
            self._create()
            return self.annotations
        raise AttributeError(name)

    def annotate(self,
            attr_name,
            annotate_as=None,
            datatype_hint=None,
            namespace_map=None,
            namespace_key=None):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        if annotate_as is None:
            annotate_as = attr_name
        if not hasattr(self, attr_name):
            raise AttributeError(attr_name)
        value = getattr(self, attr_name)
        return self.store_annotation(
                key=annotate_as,
                value=value,
                datatype_hint=datatype_hint,
                namespace_map=namespace_map,
                namespace_key=namespace_key,
                )

    def unannotate(self, attr_name):
        """
        Remove an attribute from the list of attributes to be
        persisted as an annotation.
        """
        del self.annotations.pop[attr_name]

    def clear_annotations(self):
        """
        Clears registry of annotations to be persisted.
        """
        self.annotations.clear()

    def has_annotations(self):
        """
        Returns True if there are attributes to be persisted as
        annotations.
        """
        return bool(len(self.annotations) > 0)

class Annotation(Annotated):
    "Tracks the basic information need to serialize an attribute correctly."
    def __init__(self,
            key,
            value,
            datatype_hint=None,
            namespace_map=None,
            namespace_key=None
            ):
        self.key = key
        self.value = value
        self.datatype_hint = datatype_hint
        self._namespace_map = None
        self.namespace_map = namespace_map
        self._namespace_key = None
        self.namespace_key = namespace_key

    def _get_namespace_map(self):
        if self._namespace_map is None:
            self._namespace_map = {"http://packages.python.org/DendroPy/": "dendropy"}
        return self._namespace_map
    def _set_namespace_map(self, value):
        self._namespace_map = value
    namespace_map = property(_get_namespace_map, _set_namespace_map)

    def _get_namespace_key(self):
        if self._namespace_key is None:
            self._namespace_key = "dendropy"
    def _set_namespace_key(self, value):
        self._namespace_key = value
    namespace_key = property(_get_namespace_key, _set_namespace_key)

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
