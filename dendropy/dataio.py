#! /usr/bin/env python

############################################################################
##  dataio.py
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
Data i/o.
"""

import sys
import StringIO

import dendropy
from dendropy import datasets
from dendropy import trees
from dendropy import taxa
from dendropy import nexml
from dendropy import nexus
from dendropy import phylip
from dendropy import fasta

## file formats ##
NEXUS       = 100
NEWICK      = 200
NEXML       = 300
PHYLIP      = 400
FASTA       = 600

def get_dataset(reader, file=None, string=None):
    """
    Convenience wrapper around reader.   
    """
    if file is None and string is None:
        raise Exception("File or string source must be specified.")            
    if isinstance(file, str):
        file = open(file, "r")
    elif string is not None:
        file = StringIO.StringIO(string)        
    return reader.read_dataset(file)
    
def dataset_from_nexml(file=None, string=None):
    """
    Reads a nexml file and returns a corresponding Dataset object.
    """
    get_dataset(reader=nexml.NexmlReader(), file=file, string=string)

def dataset_from_newick(file=None, string=None):
    """
    Reads a nexml file and returns a corresponding Dataset object.
    """
    get_dataset(reader=nexus.NewickReader(), file=file, string=string)

def dataset_from_nexus(file=None, string=None):
    """
    Reads a nexml file and returns a corresponding Dataset object.
    """
    get_dataset(reader=nexus.NexusReader(), file=file, string=string)    
    
    