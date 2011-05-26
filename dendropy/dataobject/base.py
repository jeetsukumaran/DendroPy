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
    def __init__(self, attr_name=None, type_hint=None):
        "__init__ creates an Annotation with optional initialization values."
        self.attr_name = attr_name  # name of attribute
        self.type_hint = type_hint  # override introspective type determination

class AnnotesDict(dict):
    "Dictionary that maps annotation terms to values."
    annotation_dict_instances = 0

    def __init__(self):
        AnnotesDict.annotation_dict_instances += 1
        self.oid = "dict" + str(AnnotesDict.annotation_dict_instances)

class DataObject(object):
    """
    Base class for elements that will be serialized.
    """
    def __init__(self):
        self.attributes = []    ## extra attributes to be written during XML serialization
        self.extensions = []    ## DOM extensions to be written during XML serialization

class Annotated(DataObject):
    """
    Base class from which all classes that need to persist object attributes
    beyond the core elements (such as id, label, etc.) will derive.
    """

    def __init__(self):
        "__init__ creates dictionary to track attributes that will be persisted."
        DataObject.__init__(self)
        self._annotations = {}

    def annotate(self, attr_name, annotate_as=None, type_hint=None):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        if annotate_as is None:
            annotate_as = attr_name
        if attr_name not in self._annotations:
            if not hasattr(self, attr_name):
                raise AttributeError(attr_name)
            else:
                annote = Annotation(attr_name=attr_name,
                                    type_hint=type_hint)
                self._annotations[annotate_as] = annote
        else:
            annote = self._annotations[annotate_as]
            annote.attr_name = annotate_as
            annote.type_hint = type_hint

    def unannotate(self, attr_name):
        """
        Remove an attribute from the list of attributes to be
        persisted as an annotation.
        """
        self._annotations.pop(attr_name, None)

    def annotations(self):
        """
        Returns the internal dictionary of annotations as a
        dictionary. If values are of Annotable class, then their
        'annotations()' method is called to populate the dictionary
        value.
        """
        annote_dict = AnnotesDict()
        for key, value in self._annotations.items():
            if hasattr(self, value.attr_name):
                attr_value = getattr(self, value.attr_name)
                if isinstance(attr_value, Annotated):
                    annote_dict[key] = attr_value.as_annotation()
                else:
                    annote_value = attr_value
                    annote_dict[key] = (annote_value, value.type_hint)
            else:
                pass # ignore deleted or non-existing attributesXS
        return annote_dict

    def as_annotation(self):
        """
        Returns a *pair of values*, self as annotation dictionary, and
        type hint if to be serialized as other than a dict. Called
        when an Annotable object is itself an annotation. Deriving
        classes should override this to return whatever information
        they want to persist when used as an annotation of another
        Annotable object.
        """
        return self.annotations(), None

    def clear_annotations(self):
        """
        Clears registry of annotations to be persisted.
        """
        self._annotations.clear()

    def has_annotations(self):
        """
        Returns True if there are attributes to be persisted as
        annotations.
        """
        return bool(len(self._annotations) > 0)

class Labelled(Annotated):
    "Provides for getting and setting of an object label."

    def __init__(self, label=None):
        """
        __init__ calls Annotated.__init__, and then, if keyword
        argument `label` is given, assigns it to self.label.
        """
        Annotated.__init__(self)
        self.label = label

class IdTagged(Labelled):
    """
    Provides infrastructure for the maintenance of a unique object id,
    ensuring that this will never be None.
    """

    instances = 0

    def normalize_id(id_str):
        """
        Given a string `id_str`, this returns a xs:NCName compliant
        version of the string: (Letter | '_' | ':')
        (NameChar)*. NameChar is given by : Letter | Digit | '.' | '-'
        | '_' | ':'
        """
        if len(id_str) > 0:
            f = id_str[0]
        else:
            f = '_' + str(id(id_str))
        if not (f.isalpha or f == '_' or f == ':'):
            id_str = '_' + id_str
        id_str = re.sub('[^\w\d\-\.]', '', id_str)
        return id_str

    normalize_id = staticmethod(normalize_id)

    def __init__(self, **kwargs):
        """
        __init__ calls Labelled.__init__, and assigns element id if
        given.
        """
        Labelled.__init__(self, label=kwargs.get("label", None))
        IdTagged.instances += 1
        if "oid" not in kwargs or kwargs["oid"] is None:
            self._oid = self._default_oid()
        else:
            self._oid = self.normalize_id(kwargs["oid"])

    def _default_oid(self):
        "Returns default oid."
        return self.__class__.__name__ + str(id(self))

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
            self._oid = self.normalize_id(oid)
        else:
            self._oid = self._default_oid()

    oid = property(_get_oid, _set_oid)

    def __str__(self):
        "String representation of the object: it's id."
        return str(self.oid)
