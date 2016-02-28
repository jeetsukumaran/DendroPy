#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
XML-parsing abstraction layer.
"""

from xml.etree import ElementTree
from dendropy.utility.textprocessing import StringIO

class XmlNamespaces(object):

    def __init__(self):
        self.namespace_prefix_map = {}
        self.prefix_namespace_map = {}

    def add_namespace(self, prefix, namespace):
        try:
            self.namespace_prefix_map[namespace].append(prefix)
        except KeyError:
            self.namespace_prefix_map[namespace] = [prefix]
        self.prefix_namespace_map[prefix] = namespace

class XmlObject(object):

    def __init__(self):
        pass

    def recast_element(self, element, subelement_factory=None):
        if element is None:
            return None
        if subelement_factory is not None:
            return subelement_factory(element)
        else:
            return self.__class__(element)

    def getiterator(self, tag, subelement_factory=None):
        for element in self._element.getiterator(tag):
            yield self.recast_element(element=element, subelement_factory=subelement_factory)

    def findall(self, tag, subelement_factory=None):
        for element in self._element.findall(tag):
            yield self.recast_element(element, subelement_factory=subelement_factory)

    def find(self, tag, subelement_factory=None):
        return self.recast_element(self._element.find(tag), subelement_factory=subelement_factory)

    def get(self, attrib_name, default=None):
        return self._element.get(attrib_name, default)

    def _get_attrib(self):
        return self._element.attrib
    attrib = property(_get_attrib)

    def _get_text(self):
        return self._element.text
    text = property(_get_text)

class XmlElement(XmlObject):
    """
    Abstraction layer around an item.
    """

    def __init__(self, element, default_namespace=None):
        XmlObject.__init__(self)
        self._element = element
        self.default_namespace = default_namespace

    def subelement_factory(self, element):
        return self.__class__(element, default_namespace=self.default_namespace)

    def format_namespace(self, namespace=None):
        if namespace:
            return "{%s}" % namespace
        if self.default_namespace:
            return "{%s}" % self.default_namespace
        return ""

    def compose_tag(self, tag, namespace=None):
        if namespace:
            ns = "{%s}" % namespace
        elif self.default_namespace:
            ns = "{%s}" % self.default_namespace
        else:
            ns = ""
        if isinstance(tag, list):
            return "/".join( ("%s%s" % (ns, i)) for i in tag )
        else:
            return "%s%s" % (ns, tag)

    def namespaced_getiterator(self, tag, namespace=None, subelement_factory=None):
        if subelement_factory is None:
            subelement_factory = self.subelement_factory
        for element in self._element.getiterator(self.compose_tag(tag, namespace)):
            yield self.recast_element(element=element, subelement_factory=subelement_factory)

    def namespaced_findall(self, tag, namespace=None, subelement_factory=None):
        if subelement_factory is None:
            subelement_factory = self.subelement_factory
        for element in self._element.findall(self.compose_tag(tag, namespace)):
            yield self.recast_element(element=element, subelement_factory=subelement_factory)

    def namespaced_find(self, tag, namespace=None, subelement_factory=None):
        if subelement_factory is None:
            subelement_factory = self.subelement_factory
        e = self._element.find(self.compose_tag(tag, namespace))
        return self.recast_element(e, subelement_factory=subelement_factory)

    def namespaced_findtext(self, tag, namespace=None):
        return self._element.findtext(self.compose_tag(tag, namespace))

class XmlDocument(XmlObject):
    """
    Abstraction layer around an XML document.
    """

    def __init__(self,
            file_obj=None,
            subelement_factory=None):
        """
        __init__ initializes a reference to the ElementTree parser, passing it
        the a file descripter object to be read and parsed or the
        ElemenTree.Element object to be used as the root element.
        """
        XmlObject.__init__(self)
        if subelement_factory is None:
            self.subelement_factory = XmlElement
        else:
            self.subelement_factory = subelement_factory
        self.namespace_registry = XmlNamespaces()
        self.root = None
        if file_obj:
            self.parse_file(file_obj)

    def parse_string(self, source):
        "Loads an XML document from an XML string, source."
        s = StringIO(source)
        return self.parse_file(source)

    def parse_file(self, source):
        """
        Loads an XML document from source, which can either be a
        filepath string or a file object.
        Custom parsing to make sure namespaces are saved.
        """
        events = "start", "start-ns", "end-ns"
        root = None
        ns_map = []
        for event, elem in ElementTree.iterparse(source, events):
            if event == "start-ns":
                ns_map.append(elem)
            elif event == "start":
                if root is None:
                    root = elem
        # self.doctree = ElementTree.ElementTree(element=root)
        self.root = self.subelement_factory(root)
        for prefix, namespace in ns_map:
            self.namespace_registry.add_namespace(prefix=prefix, namespace=namespace)

