#! /usr/bin/env python

############################################################################
##  xmlparser.py
##
##  Part of the DendroPy phylogenetic computation library.
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
This module abstracts the parsing of XML.
"""

import sys
if sys.version_info[0] >= 2 and sys.version_info[1] >= 5:
    from xml.etree import ElementTree
else:
    try:
        from elementtree import ElementTree
    except ImportError:
        print >> sys.stderr, \
              'XML parsing requires the ElementTrees package '\
              '(http://effbot.org/zone/element.htm) installed ' \
              'if using Python older than Version 2.5.0.'
        raise ImportError("ElementTrees package not available")

from dendropy import utils

class xml_document(object):
    """
    ElementTree requires that the complete XML be loaded in memory
    before working with it, which may be discouraging when dealing
    with large files. By abstracting it this way, if we want to
    implement a more efficient approach later, we will avoid the need
    for messing with other sections of code.
    """

    def __init__(self, element=None, file=None):
        """
        Initializes reference to the ElementTree parser, passing it
        the a file descripter object to be read and parsed or the
        ElemenTree.Element object to be used as the root element.
        """
        self.etree = ElementTree.ElementTree(element=element, file=file)

    def parse_string(self, source):
        """
        Loads an XML document from an XML string, source.
        """
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
        return utils.RecastingIterator(self.etree.getroot().getiterator(tag), \
                                       XmlElement)
        
class XmlElement(object):
    """
    Generic XML element. May contain child elements. At present a
    simple wrapper around the ElementTree.Element class, and
    implementing just the components of the ElementTree.Element
    interface that are needed for DendroPy.
    """

    def __init__(self, element):
        """
        Initializes basic structure of object. `element` is an
        ElementTree.Element object.
        """
        self.etree_element = element

    def getiterator(self, tag):
        """
        Returns an iterator over child elements with tags that match `tag`.
        """
        return utils.RecastingIterator(self.etree_element.getiterator(tag), \
                                       XmlElement)

    def get(self, key, default=None):
        """
        Returns the attribute of this element with matching key, or
        substituting default if not found.
        """
        return self.etree_element.get(key, default)

    def find(self, path):
        """
        Finds all matching subelements, by tag name or path.
        """
        return self.etree_element.find(path)
