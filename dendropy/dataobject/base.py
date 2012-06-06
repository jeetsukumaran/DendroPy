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

class Annotation(object):
    "Tracks the basic information need to serialize an attribute correctly."
    def __init__(self,
            annotate_as,
            value,
            datatype_hint=None,
            namespace=None,
            namespace_prefix=None
            ):
        self.annotate_as = annotate_as
        self.value = value
        self.datatype_hint = datatype_hint
        self._namespace = None
        self.namespace = namespace
        self._namespace_prefix = None
        self.namespace_prefix = namespace_prefix

    def _get_namespace(self):
        if self._namespace is None:
            self._namespace = "http://packages.python.org/DendroPy/"
        return self._namespace

    def _set_namespace(self, value):
        self._namespace = value

    namespace = property(_get_namespace, _set_namespace)

    def _get_namespace_prefix(self):
        if self._namespace_prefix is None:
            self._namespace_prefix = "dendropy"

    def _set_namespace_prefix(self, value):
        self._namespace_prefix = value

    namespace_prefix = property(_get_namespace_prefix, _set_namespace_prefix)


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
        self.annotations = {}

    def store_annotation(self,
            annotate_as,
            value,
            datatype_hint=None,
            namespace=None,
            namespace_prefix=None):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        if not hasattr(self, "_annotations"):
            self._create()
        annote = Annotation(
                annotate_as=annotate_as,
                value=value,
                datatype_hint=datatype_hint,
                namespace=namespace,
                namespace_prefix=namespace_prefix,
                )
        self.annotations[annotate_as] = annote
        return annote

    def annotate(self,
            attr_name,
            annotate_as=None,
            datatype_hint=None,
            namespace=None,
            namespace_prefix=None):
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
                annotate_as=annotate_as,
                value=value,
                datatype_hint=datatype_hint,
                namespace=namespace,
                namespace_prefix=namespace_prefix,
                )

    def unannotate(self, attr_name):
        """
        Remove an attribute from the list of attributes to be
        persisted as an annotation.
        """
        if not hasattr(self, "_annotations"):
            self._create()
        del self._annotations.pop[attr_name]

    def clear_annotations(self):
        """
        Clears registry of annotations to be persisted.
        """
        if not hasattr(self, "_annotations"):
            self._create()
        self._annotations.clear()

    def has_annotations(self):
        """
        Returns True if there are attributes to be persisted as
        annotations.
        """
        if not hasattr(self, "annotations"):
            self._create()
        return bool(len(self.annotations) > 0)

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
