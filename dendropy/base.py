#! /usr/bin/env python

############################################################################
##  base.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Infrastructure for object serialization and management. 
"""

import re

class Annotation(object):
    "Tracks the basic information need to serialize an attribute correctly."
    def __init__(self, attr_name=None, type_hint=None):
        "Instantiates with optional initialization values."
        self.attr_name = attr_name  # name of attribute
        self.type_hint = type_hint  # override introspective type determination

class AnnotesDict(dict):
    annotation_dict_instances = 0
    
    def __init__(self):
        AnnotesDict.annotation_dict_instances += 1
        self.oid = "dict" + str(AnnotesDict.annotation_dict_instances)

class Element(object):
    """
    Base class for elements that will be serialized.
    """
    def __init__(self):
        self.attributes = []    ## extra attributes to be written during XML serialization
        self.extensions = []    ## DOM extensions to be written during XML serialization
        
class Annotated(Element):
    """
    Base class from which all classes that need to persist object attributes
    beyond the core elements (such as id, label, etc.) will derive.
    """
    
    def __init__(self):
        "Creates dictionary to track attributes that will be persisted."
        Element.__init__(self)
        self.__annotations = {}

    def annotate(self, attr_name, annotate_as=None, type_hint=None):
        """
        Add an attribute to the list of attributes that need to be
        persisted as an annotation.
        """
        #@ Jeet, this next line used to be
        #annotate_as = annotate_as is None and attr_name or attr_name
        # I assume that was a bug (although it was not causing tests to fail).
        annotate_as = annotate_as is None and attr_name or annotate_as  
        if attr_name not in self.__annotations:
            if not hasattr(self, attr_name):
                raise AttributeError(attr_name)
            else:
                annote = Annotation(attr_name=attr_name,
                                    type_hint=type_hint)
                self.__annotations[annotate_as] = annote
        else:
            annote = self.__annotations[annotate_as]
            annote.attr_name = annotate_as
            annote.type_hint = type_hint
            
    def unannotate(self, attr_name):
        """
        Remove an attribute from the list of attributes to be
        persisted as an annotation.
        """
        self.__annotations.pop(attr_name, None)

    def annotations(self):
        """
        Returns the internal dictionary of annotations as a
        dictionary. If values are of Annotable class, then their
        'annotations()' method is called to populate the dictionary
        value.
        """        
        annote_dict = AnnotesDict()
        for key, value in self.__annotations.items():
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

    def has_annotations(self):
        """
        Returns True if there are attributes to be persisted as
        annotations.
        """
        return bool(len(self.__annotations) > 0)

class Labelled(Annotated):
    "Provides for getting and setting of an object label."

    def __init__(self, label=None):
        """
        Initializes by calling base class, and then, if keyword
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
        f = id_str[0]
        if not (f.isalpha or f == '_' or f == ':'):
            id_str = '_' + id_str
        id_str = re.sub('[^\w\d\-\.]', '', id_str)
        return id_str

    normalize_id = staticmethod(normalize_id)

    def __init__(self, oid=None, label=None):
        """
        Initializes by calling base classes, and assigns element id if
        given.
        """
        Labelled.__init__(self, label=label)
        IdTagged.instances += 1
        if oid is None:
            self.__oid = self._default_oid()
        else:
            self.__oid = self.normalize_id(oid)
            
    def _default_oid(self):
        "Returns default oid."
        if self.label is None:
            prefix = self.__class__.__name__
            return "_" + prefix + str(IdTagged.instances)            
        else:
            return self.normalize_id(self.label)        
                
    def _get_oid(self):
        "Returns id."
        return self.__oid

    def _set_oid(self, oid):
        """
        Sets oid to oid if oid is not None (normalized to conform
        to xs:NCName specs if neccessary), otherwise sets to some
        other oidue.
        """
        if oid is not None:
            self.__oid = self.normalize_id(oid)
        else:
            self.__oid = self._default_oid()

    oid = property(_get_oid, _set_oid)
