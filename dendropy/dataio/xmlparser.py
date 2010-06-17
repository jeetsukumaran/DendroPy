#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

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

nexml_namespace = "http://www.nexml.org/1.0"
diagnosed_tags = []


def diagnose_namespace(tag, namespace):
    if tag not in diagnosed_tags:
        diagnosed_tags.append(tag)
        sys.stderr.write("% 20s\t%s\n" % (tag, namespace))


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
            print "decorated_tag = ", decorated_tag
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
    if i:
        diagnose_namespace(tag, "no namespace decoration")
    elif namespace_list:
        d = {'tag' : tag}
        for n in namespace_list:
            d['ns'] = n
            i = meth("{%(ns)s}$(tag)s" % d)
            if i:
                diagnose_namespace(tag, "decorated with namespace %(ns)s" % d)
                break
    if not i:
        diagnose_namespace(tag, "NOT FOUND")
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
