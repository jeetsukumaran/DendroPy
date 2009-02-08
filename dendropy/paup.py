#! /usr/bin/env python

############################################################################
##  paup.py
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

"""PAUP* wrapping/parsing"""


import os
import subprocess
import tempfile
import re

from dendropy import taxa

if "PAUP_PATH" in os.environ:
    PAUP_PATH = os.environ["PAUP_PATH"]
else:
    PAUP_PATH = "paup"

class Paup(object):
    """ Wrapper around PAUP* """
    
    def __init__(self, paup_path=None):
        if paup_path is None:
            self.paup_path = PAUP_PATH
        else:
            self.paup_path = paup_path
            
    def run(self, commands):
        """ executes given list of commands in PAUP*, 
        return results of stdout and stderr """
        commands = "\n".join(commands) + "\n"
        paup_run = subprocess.Popen(['%s -n' % self.paup_path],
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
        stdout, stderr = paup_run.communicate(commands)
        return stdout.split("\n")
        
    def parse_taxa_block(self, lines):
        """
        Given PAUP* output that includes a taxon listing as produced by
        `compose_list_taxa`, this parses out and returns a taxon block.
        """
        taxlabels = []
        taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
        idx = 0
        for line in lines:
            idx += 1
            if line == "TAXON LIST BEGIN":
                break                  
        for line in lines[idx:]:
            if line == "TAXON LIST END":
                break
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                taxlabels.append(ti_match.group(2).strip())                
        taxa_block = taxa.TaxaBlock() 
        for taxlabel in taxlabels:
            taxa_block.add_taxon(label=taxlabel.replace(' ', '_'))
        return taxa_block            
        
    def compose_list_taxa(self):
        """ 
        Given a data file in memory, this gets PAUP* to print a list of 
        taxa that can be used to build a TaxaBlock later.
        """
        return ["[!TAXON LIST BEGIN]\ntstatus / full;\n[!TAXON LIST END]\n"]
    
    def compose_load_trees(self,
                           tree_filepaths,
                           taxa_filepath=None, # for taxa block; leave None if taxa in treefile                                                     
                           burnin=0,
                           mode=7, # keep trees in memory, specify 3 to clear
                           reset=True):
        """
        Composes commands to load a set of trees into PAUP*, with the specified 
        number of burnin dropped.
        """
        if isinstance(tree_filepaths, str):
            raise Exception("expecting list of filepaths, not string")
        gettree_template = 'gett file= %%s storebrlens=yes warntree=no unrooted=yes from=%d mode=%%d;' % (burnin+1)
        paup_template = []
        paup_template.append("set warnreset=no; set increase=auto; set warnroot=no;")
        if taxa_filepath is not None:
            paup_template.append('execute %s;' % taxa_filepath)
            paup_template.append(gettree_template % (tree_filepaths[0], mode))
        else:
            if reset:
                paup_template.append('execute %s;' % tree_filepaths[0])
            else:
                paup_template.append(gettree_template % (tree_filepaths[0], mode))
        for tree_filepath in tree_filepaths[1:]:
            paup_template.append(gettree_template % (tree_filepath, 7))
                        
        return paup_template                                

if __name__ == "__main__":
    paup = Paup()
    treefiles = ["/Users/jeet/Documents/Projects/Phyloinformatics/DendroPy/dendropy/dendropy/tests/data/feb032009.tre"] #, "/Users/jeet/Documents/Projects/Phyloinformatics/DendroPy/dendropy/dendropy/tests/data/feb032009.tre"]
    stdout, stderr = paup.run(paup.compose_load_trees(treefiles) + paup.compose_list_taxa())
    print stdout