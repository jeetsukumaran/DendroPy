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
XML-parsing abstraction layer.
"""

import sys
if sys.version_info[0] >= 2 and sys.version_info[1] >= 5:
    from xml.etree import ElementTree
else:
    try:
        from elementtree import ElementTree
    except ImportError:
        try:
            from dendropy.utility import ElementTree
        except ImportError:
            sys.stderr.write("""\

    ###############################################################################
    Failed to import the XML parsing module.
    If you are trying to install DendroPy, then the installation failed because
    Python versions of less than 2.5.0 require the ElementTrees package installed.

    Please either use a newer version of Python 2 (e.g., Python 2.5, 2.6 or 2.7) to
    install DendroPy, or install the ElementTrees package from:

        http://effbot.org/zone/element.htm
        http://effbot.org/downloads/#elementtree
    ###############################################################################

    """)
            sys.exit(1)

from dendropy.utility import containers

diagnosed_tags = []


def diagnose_namespace(tag, namespace):
    if tag not in diagnosed_tags:
        diagnosed_tags.append(tag)
#        sys.stderr.write("% 20s\t%s\n" % (tag, namespace))


def _getiterator(etree, tag, namespace_list=()):
    """
    Returns an iterator over all top-level elements from the root element
    that have the matching tag.
    """
    i = etree.getiterator(tag)
    if i:
        diagnose_namespace(tag, "no namespace decoration")
    elif namespace_list:
        d = {'tag' : tag}
        for n in namespace_list:
            d['ns'] = n
            decorated_tag = "{%(ns)s}%(tag)s" % d
            #print "decorated_tag = ", decorated_tag
            i = etree.getiterator(decorated_tag)
            if i:
                diagnose_namespace(tag, "decorated with namespace %(ns)s" % d)
                break
    if not i:
        diagnose_namespace(tag, "NOT FOUND")
    recasting = lambda x: XmlElement(x, namespace_list=namespace_list)
    return containers.RecastingIterator(i, recasting)


def _invoke_method_for_namespaces(meth, tag, namespace_list=()):
    i = meth(tag)
    if i is not None:
        diagnose_namespace(tag, "no namespace decoration")
    elif namespace_list:
        d = {'tag' : tag}
        for n in namespace_list:
            d['ns'] = n
            decorated_tag = "{%(ns)s}%(tag)s" % d
            #print "decorated_tag = ", decorated_tag
            i = meth(decorated_tag)
            if i is not None:
                diagnose_namespace(tag, "decorated with namespace %(ns)s" % d)
                break
    if i is None:
        diagnose_namespace(tag, "NOT FOUND")
        return i
    if not isinstance(i, str) and not isinstance(i, unicode):
        return XmlElement(i, namespace_list=namespace_list)
    return i

class xml_document(object):
    """
    ElementTree requires that the complete XML be loaded in memory
    before working with it, which may be discouraging when dealing
    with large files. By abstracting it this way, if we want to
    implement a more efficient approach later, we will avoid the need
    for messing with other sections of code.
    """

    def __init__(self, element=None, file_obj=None, namespace_list=()):
        """
        __init__ initializes a reference to the ElementTree parser, passing it
        the a file descripter object to be read and parsed or the
        ElemenTree.Element object to be used as the root element.
        """
        self.namespace_list = list(namespace_list)
        self.etree = ElementTree.ElementTree(element=element, file=file_obj)

    def parse_string(self, source):
        "Loads an XML document from an XML string, source."
        root = ElementTree.fromstring(source)
        self.etree = ElementTree(element=root)

    def parse_file(self, source):
        """
        Loads an XML document from source, which can either be a
        filepath string or a file object.
        """
        root = ElementTree.parse(source=source)
        self.etree = ElementTree(element=root)

    def getiterator(self, tag):
        """
        Returns an iterator over all top-level elements from the root element
        that have the matching tag.
        """
        return _getiterator(self.etree.getroot(), tag, self.namespace_list)

class XmlElement(object):
    """
    Generic XML element. May contain child elements. At present a
    simple wrapper around the ElementTree.Element class, and
    implementing just the components of the ElementTree.Element
    interface that are needed for DendroPy.
    """

    def __init__(self, element, namespace_list=()):
        """
        __init__ initializes a basic structure of object. `element` is an
        ElementTree.Element object.
        """
        self.namespace_list = namespace_list
        self.etree_element = element

    def getiterator(self, tag):
        "Returns an iterator over child elements with tags that match `tag`."
        return _getiterator(self.etree_element, tag, self.namespace_list)

    def get(self, key, default=None):
        """
        Returns the attribute of this element with matching key, or
        substituting default if not found.
        """
        i = _invoke_method_for_namespaces(self.etree_element.get, key, self.namespace_list)
        if not i:
            return default
        return i

    def findtext(self, text):
        "Finds free text contained in element"
        return _invoke_method_for_namespaces(self.etree_element.findtext, text, self.namespace_list)

    def find(self, path):
        "Finds all matching subelements, by tag name or path."
        return _invoke_method_for_namespaces(self.etree_element.find, path, self.namespace_list)
