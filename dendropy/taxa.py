#! /usr/bin/env python

############################################################################
##  taxa.py
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
This module provides classes and methods for managing taxa.
"""
import sys
from dendropy import base


_HEX_CODE_TO_BINARY = {'a': '1010', 'c': '1100', 'b': '1011', 'e': '1110', 'd': '1101', 'f': '1111', '1': '0001', '0': '0000', '3': '0011', '2': '0010', '5': '0101', '4': '0100', '7': '0111', '6': '0110', '9': '1001', '8': '1000'}

def int_to_bitstring(n):
    global _HEX_CODE_TO_BINARY
    assert n >= 0
    if n == 0:
        return "0"
    h = hex(n)
    b = "".join([_HEX_CODE_TO_BINARY[i] for i in h[2:]])
    return b[b.find("1"):] 

class TaxonLinked(base.IdTagged):
    """
    Provides infrastructure for maintaining link/reference to a Taxon
    object.
    """
    
    def __init__(self, oid=None, label=None, taxon=None):
        """
        Initializes by calling base class.
        """
        base.IdTagged.__init__(self, oid=oid, label=label)
        self.__taxon = taxon

    def _get_taxon(self):
        """
        Returns taxon associated with this object.
        """
        return self.__taxon

    def _set_taxon(self, taxon):
        """
        If `taxon` is a Taxon object, then it is assigned directly. If
        `taxon` is a string, then it is assumed to be a label, and a
        new taxon object is constructed based on it and assigned (the
        new taxon object will have the string given by `taxon` as id
        and label, though it will be modified as neccessary to make it
        xs:NCName compliant for the id.
        """
        if taxon is None:
            self.__taxon = None
        elif isinstance(taxon, Taxon):
            self.__taxon = taxon
        else:
            taxon_obj = Taxon()
            taxon_obj.label = taxon
            taxon_obj.oid = taxon
            self.__taxon = taxon_obj

    taxon = property(_get_taxon, _set_taxon)

class TaxaLinked(base.IdTagged):
    """
    Provides infrastructure for the maintenance of references to taxa
    blocks.
    """

    def __init__(self, oid=None, label=None, taxa_block=None):
        """
        Initializes by calling base class.
        """
        base.IdTagged.__init__(self, oid=oid, label=label)
        self.__taxa_block = taxa_block

    def _get_taxa_block(self):
        """
        Returns taxon block associated with this object. If none is
        given, then it builds one.
        """
        if self.__taxa_block is None:
            self.__taxa_block= TaxaBlock()
        return self.__taxa_block

    def _set_taxa_block(self, taxa_block):
        """
        Sets the taxon block for this object.
        """
        self.__taxa_block = taxa_block

    taxa_block = property(_get_taxa_block, _set_taxa_block)
            
class TaxaBlock(list, base.IdTagged):
    """
    Taxon manager.
    """

    def __init__(self, *args, **kwargs):
        """
        Inits. Handles keyword arguments: `oid` and `label`.
        """
        list.__init__(self, *args)
        base.IdTagged.__init__(self, *args, **kwargs)

    def __str__(self):
        """
        String representation of self.
        """
        header = []
        if self.oid:
            header.append("%s" % str(self.oid))
        if self.label:
            header.append("(\"%s\")" % self.label)
        taxlist = []
        for taxon in self:
            taxlist.append(str(taxon))
        return ' '.join(header) + ' : [' + ', '.join(taxlist) + ']' 

    def get_taxon(self, oid=None, label=None, update=False):
        """
        Retrieves taxon object with given id OR label (if both are
        given, the first match found is returned). If taxon does not
        exist and update is False, an exception is raised. If taxon
        does not exist and update is True, then a new taxon is
        created, added, and returned.
        """
        if not oid and not label:
            raise Exception("Need to specify Element ID or Label.")
        for taxon in self:
            if taxon.oid == oid \
               or taxon.label == label:
                return taxon
        if not update:
            raise Exception("Taxon not found: %s/%s" % (oid, label))
        else:
            taxon = Taxon(oid=oid, label=label)
            self.append(taxon)
            return taxon
            
    def add_taxon(self, oid=None, label=None):
        """
        Convenience function that wraps `get_taxon`.
        """
        self.get_taxon(oid=oid, label=label, update=True)
        
    def clear(self):
        """
        Removes all taxa from this block.
        """
        for t in self:
            self.remove(t)
            
    def labels(self):
        """
        Convenience method to return all taxa labels.
        """
        return [str(taxon.label) for taxon in self]        
        
    def complement_split_bitmask(self, split):
        """
        Returns complement of the split bitmask.
        """
        return split ^ self.all_taxa_bitmask()
            
    def all_taxa_bitmask(self):
        """
        Returns mask of all taxa.
        """
        return pow(2, len(self)) - 1
        
    def taxon_bitmask(self, taxon):
        """
        Returns unique bitmask of given taxon. Will raise index error if
        taxon does not exist.
        """
        try:
            return pow(2, self.index(taxon))
        except ValueError:
            raise ValueError("Taxon with ID '%s' and label '%s' not found" 
                             % (str(taxon.oid), str(taxon.label)))        
                             
    def split_bitmask_string(self, split_bitmask):
        """
        Returns bitstring representation of split_bitmask.
        """
        return "%s" % int_to_bitstring(split_bitmask).rjust(len(self), "0")
                                
class Taxon(base.IdTagged):
    """
    A taxon associated with a sequence or a node on a tree.
    """
    
    def cmp(taxon1, taxon2):
        """
        Compares taxon1 and taxon2 based on label.
        """
        return cmp(str(taxon1.label), str(taxon2.label))
    
    cmp = staticmethod(cmp)

    def __init__(self, oid=None, label=None): 
        """
        Initializes by calling base class.
        """
        base.IdTagged.__init__(self, oid=oid, label=label)

    def __str__(self):
        """
        String representation of self = taxon name.
        """
        return str(self.label)
