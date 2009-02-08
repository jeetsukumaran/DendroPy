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

"""
Various calls to PAUP* to calculate stuff.
"""

import os
import sys
import subprocess
import tempfile
import re

import unittest
import dendropy.tests
from dendropy import get_logger
_LOG = get_logger("PAUPWrapper")

from dendropy import datasets
from dendropy import dataio
from dendropy import taxa
from dendropy import splits

###############################################################################
## PAUP* WRAPPER
###############################################################################

if "PAUP_PATH" in os.environ:
    PAUP_PATH = os.environ["PAUP_PATH"]
else:
    PAUP_PATH = "paup"

class Paup(object):
    """ Wrapper around PAUP* """
    
    def paup_group_to_mask(group_string):
        """
        This converts a PAUP* group representation (i.e. a string of askterisks
        and periods, where the asterisks denote the taxon index counting from
        left to right) to a mask representation:
            - a clade mask, where 1's represent descendents of the split/edge 
              (with taxon index counting from right to left, i.e., first taxon
              is right-most bit)
            - a split mask, an unrooted normalized version of the above, where 
              if the right most bit is not 1 the clade mask is complemented 
              (and not changed otherwise).
        """
        group_string = group_string[::-1] # flip to get correct orientation
        bit_string = group_string.replace("*", "1").replace(".", "0")
        clade_mask = int(bit_string, 2)
        if not (clade_mask & 1):
            split_mask = int(group_string.replace("*", "0").replace(".", "1"), 2)
        else:
            split_mask = clade_mask[:]
        return clade_mask, split_mask
        
    paup_group_to_mask = staticmethod(paup_group_to_mask)        
    
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
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        stdout, stderr = paup_run.communicate(commands)
        results = stdout.split("\n")
        if stderr:
            sys.stderr.write(stderr)
            sys.exit(1)
        # check results:
#         for idx, line in enumerate(results):
#             try:
#                 assert not re.match("^Error\(.*", line)
#             except AssertionError, e:
#                 sys.stderr.write("Error in PAUP run:\n")
#                 # write out PAUP error message with some context
#                 i = 0
#                 while (idx+i) < len(results) and i < 3:
#                     sys.stderr.write(results[idx+i] + "\n")
#                     i += 1
        return results            
        
    def parse_taxa_block(self, paup_output):
        """
        Given PAUP* output that includes a taxon listing as produced by
        `compose_list_taxa`, this parses out and returns a taxon block.
        """
        taxlabels = []
        taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
        idx = 0
        for line in paup_output:
            idx += 1
            if line == "TAXON LIST BEGIN":
                break                  
        for line in paup_output[idx:]:
            if line == "TAXON LIST END":
                break
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                taxlabels.append(ti_match.group(2).strip())                
        taxa_block = taxa.TaxaBlock() 
        for taxlabel in taxlabels:
            taxa_block.add_taxon(label=taxlabel.replace(' ', '_'))
        return taxa_block           
        
    def parse_group_freqs(self, paup_output):
        """
        Given PAUP* output that includes a split counting procedure,
        this collects the splits and returns a dictionary of group strings and 
        their frequencies
        """        
        bipartitions = []
        bipartition_freqs = {}
        bipartition_counts = {}
        bipartition_pattern = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')       
        idx = 0
        for line in paup_output:
            idx += 1
            if line == "SPLITS COUNT BEGIN":
                break                  
        for line in paup_output[idx:]:
            if line == "SPLITS COUNT END":
                break
            bp_match = bipartition_pattern.match(line)
            if bp_match:
                bipartitions.append(bp_match.group(1))
                bipartition_counts[bp_match.group(1)] = int(bp_match.group(2))
                bipartition_freqs[bp_match.group(1)] = float(bp_match.group(3))
        return bipartition_counts, bipartition_freqs
        
    def split_counts_from_group_freqs(self, 
        bipartition_counts,
        taxa_block):
        """
        Returns a populated NormalizingBitmaskDict object based on the given
        bipartition info, with keys = normalized splits and values = tuple 
        (count, percentage).
        """
        #assert len(bipartition_freqs[0]) == len(taxa_block)
        pass
        
    def compose_list_taxa(self):
        """ 
        Given a data file in memory, this gets PAUP* to print a list of 
        taxa that can be used to build a TaxaBlock later.
        """
        return ["[!TAXON LIST BEGIN]\ntstatus / full;\n[!TAXON LIST END]\n"]
        
    def compose_count_splits(self, majrule_filepath=None, majrule_freq=0.5):
        """
        Given trees in memory, this composes a command to count the split
        frequencies across the trees as well as a save the majority-rule
        consensus tree if a path is given.
        """
        percent = 100 * majrule_freq
        if majrule_filepath is not None:
            treefile = " treefile=%s replace=yes "
        else:
            treefile = ""
        paup_template = []
        paup_template.extend(self.compose_list_taxa())
        paup_template.append("[!SPLITS COUNT BEGIN]")        
        paup_template.append("contree / strict=no %s showtree=no grpfreq=yes majrule=yes percent=%d" % (treefile, percent));
        paup_template.append("[!SPLITS COUNT END]")
        return paup_template
    
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


###############################################################################
## OLD STUFF
###############################################################################

def bipartitions(data_filepath,
                 tree_filepath,
                 min_clade_freq=0.5,
                 burnin=0,
                 paup_path=PAUP_PATH):
    """
    Given a set of trees (and data file), this uses PAUP*'s contree
    command to calculate the splits (bipartitions) on the trees, as well
    as their counts and relative percentages. Returned is:

        - list of taxon labels, in order of the index assigned to them by PAUP
        - list of bipartition strings in PAUP*'s notation (e.g., "...**.*.*")
        - a dictionary with the bipartition string as a key and the count of the
          bipartition occurrence in the trees examined as values
        - a dictionary with the bipartition string as a key and the
          percentage of trees with the bipartition occurence as values.
    """

    paup_args = {
        'data_filepath': data_filepath,
        'tree_filepath': tree_filepath,
        'percent': min_clade_freq * 100,
        'burnin': burnin+1,
    }
    paup_template = """\
    set warnreset=no;
    set increase=auto;
    exe %(data_filepath)s;
    gett file=%(tree_filepath)s storebrlens=yes warntree=no unrooted=yes;
    tstatus / full;
    contree %(burnin)d-. / strict=no showtree=no grpfreq=yes majrule=yes percent=%(percent)d;
"""
    paup_run = subprocess.Popen(['%s -n' % paup_path],
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    stdout, stderr = paup_run.communicate(paup_template % paup_args)
    paup_output = stdout.split('\n')
    tax_labels = []
    bipartitions = []
    bipartition_freqs = {}
    bipartition_counts = {}
    bipartition_pattern = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')
    taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
    for line in paup_output:
        bp_match = bipartition_pattern.match(line)
        if bp_match:
            bipartitions.append(bp_match.group(1))
            bipartition_counts[bp_match.group(1)] = int(bp_match.group(2))
            bipartition_freqs[bp_match.group(1)] = float(bp_match.group(3))
        else:
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                tax_labels.append(ti_match.group(2).strip())                
    return tax_labels, bipartitions, bipartition_counts, bipartition_freqs


def estimate_char_model(tree_model,
                        char_block,
                        num_states=6,
                        unequal_base_freqs=True,
                        gamma_rates=True,
                        prop_invar=True,
                        paup_path='paup'):
    """
    Returns likelihood score as well as estimates of rates, kappa, 
    base_frequencies, alpha, prop_invar, etc. (as dictionary).
    """
    tf = tempfile.NamedTemporaryFile()
    dataio.store_trees([tree_model], format='nexus', dest=tf)
    tf.flush()
    df = tempfile.NamedTemporaryFile()
    dataio.store_chars(char_block=char_block, format='nexus', dest=df)
    df.flush()    
    paup_args = {
        'datafile' : df.name,
        'treefile' : tf.name,
        'nst': num_states,
        'basefreq' : 'estimate' if unequal_base_freqs else 'equal',
        'rates' : 'gamma' if gamma_rates else 'equal',
        'pinvar' : 'estimate' if prop_invar else '0',
    }
    paup_template = """\
    set warnreset=no;
    exe %(datafile)s;
    gettrees file=%(treefile)s storebrlens=yes;
    lset tratio=estimate rmatrix=estimate nst=%(nst)s basefreq=%(basefreq)s rates=%(rates)s shape=estimate pinvar=%(pinvar)s userbrlens=yes;
    lscore 1 / userbrlens=yes;
""" 
    paup_run = subprocess.Popen(['%s -n' % paup_path],
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    stdout, stderr = paup_run.communicate(paup_template % paup_args)

    patterns = {
        'likelihood' : re.compile('-ln L\s+([\d\.]+)'),
        'AC' : re.compile('  AC\s+([\d\.]+)'),
        'AG' : re.compile('  AG\s+([\d\.]+)'),
        'AT' : re.compile('  AT\s+([\d\.]+)'),
        'CG' : re.compile('  CG\s+([\d\.]+)'),
        'CT' : re.compile('  CT\s+([\d\.]+)'),
        'GT' : re.compile('  GT\s+([\d\.]+)'),
        'kappa': re.compile('  kappa\s+([\d\.]+)'),
        'prop_invar' : re.compile('P_inv\s+([\d\.]+)'),
        'alpha' : re.compile('Shape\s+([\S]+)'),
    
    }

    results = {}
    for value_name in patterns:
        results[value_name] = None
    for line in stdout.split('\n'):
        for value_name in patterns:
            m = patterns[value_name].match(line)
            if m:
                results[value_name] = m.group(1)
                
    for value_name in results:
        if value_name == 'likelihood':
            results[value_name] = -1 * float(results[value_name])
        elif results[value_name] is not None:
            try:
                results[value_name] = float(results[value_name])
            except:
                pass
                
    return results

###############################################################################
## TEST SUITE (for new stuff)
###############################################################################
    
class PaupWrapperDumbTests(unittest.TestCase):
    """
    Checks basic running of PAUP commands, and correct parsing/extraction of 
    output. Does not do (or check) higher-level processing of output.1
    """

    def check_taxa_block(self, filename, taxlabels):
        """Loads a taxa block from `filename`, make sure taxa returned match
        `taxlabels`"""
        p = Paup()
        commands = []             
        commands.extend(p.compose_load_trees([dendropy.tests.data_source_path(filename)]))
        commands.extend(p.compose_list_taxa())
        results = p.run(commands) 
        taxa_block = p.parse_taxa_block(results)
        assert len(taxa_block) == len(taxlabels)
        for i, t in enumerate(taxa_block):
            assert t.label == taxlabels[i]
                        
    def check_group_freqs(self, treefile, group_freqs):
        """Calculates group frequencies from `filename`, make sure that 
        frequencies match `group_freqs` (given as dictionary of PAUP* group
        strings and their frequencies for the file)."""        
        p = Paup()
        commands = []             
        commands.extend(p.compose_load_trees(tree_filepaths=[dendropy.tests.data_source_path(treefile)]))
        commands.extend(p.compose_count_splits())
        results = p.run(commands)
        bipartition_counts, bipartition_freqs = p.parse_group_freqs(results)
                       
    def testTaxaBlock(self):
        test_cases = (
            ("feb032009.tre",("T01", "T02", "T03", "T04", "T05", "T06",
            "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14",
            "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22",
            "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
            "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38",
            "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46",
            "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54",
            "T55", "T56", "T57", "T58", "T59")),
            ("primates.chars.nexus", ("Lemur_catta", "Homo_sapiens",
            "Pan", "Gorilla", "Pongo", "Hylobates", "Macaca_fuscata",
            "Macaca_mulatta", "Macaca_fascicularis", "Macaca_sylvanus",
            "Saimiri_sciureus", "Tarsius_syrichta", ))
        )        
        for i in test_cases:
            self.check_taxa_block(i[0], i[1])
            
    def testGroupFreqs(self):
        test_cases = (
            ("anolis.mcmct.trees.nexus",
                {
                "......*...*..................." : (1000, 99.9), "...*...................*......" : (999, 99.8), "...................*.....*...." : (999, 99.8), "......................*.*....." : (998, 99.7), "...........*......*..........." : (994, 99.3), "............**..*..........*.." : (994, 99.3), ".*.....*......................" : (993, 99.2), "............**................" : (992, 99.1), ".***********..**.*****.*.**.**" : (991, 99.0), ".....*...*...................." : (991, 99.0), "........*......*..........*..." : (990, 98.9), "..*.........................*." : (990, 98.9), "...*.............*.....*......" : (987, 98.6), "......*.*.*....*..........*..." : (984, 98.3), "..............*......*........" : (970, 96.9), "..**.............*.....*....*." : (968, 96.7), ".***********..**.**********.**" : (961, 96.0), ".*..*..*...*......*.*........." : (894, 89.3), "........*......*.............." : (881, 88.0), "............**.............*.." : (872, 87.1), ".*..**.*.*.*..*...*.**........" : (867, 86.6),
                ".*..********..**..*.**....*..." : (707, 70.6), ".*..*..*......................" : (681, 68.0), ".*..*..*...*..*...*.**........" : (506, 50.5), ".*..*..*............*........." : (463, 46.3), "..**.............*.*...*.*..**" : (397, 39.7), "..**.............*.*...*.*..*." : (322, 32.2), ".*..********..**..*.**....*..*" : (299, 29.9), "...................*.....*...*" : (298, 29.8), ".*..*..*...*......*..........." : (295, 29.5), ".....*...*....*......*........" : (255, 25.5), "..**.............*.....*....**" : (235, 23.5), "....*...............*........." : (142, 14.2), "......*.*.*....*..........*..*" : (140, 14.0), ".*..**.*.*.*......*.*........." : (139, 13.9), ".*.....*...*......*..........." : (137, 13.7), ".***********..**.*****.*.**.*." : (128, 12.8), "...........*......*.*........." : (125, 12.5), ".*..**.*.*.*..*...*.**.......*" : (117, 11.7), "........*.................*..." : (110, 11.0), ".*..*..*...*......*.*........*" : (110, 11.0), "............**..*............." : (104, 10.4),
                "....................*........*" : (84, 8.4), ".***********..**.**.**.*..*.*." : (79, 7.9), ".***********..**.**.**.*..*.**" : (72, 7.2), ".*..********..**..****...**..." : (70, 7.0), ".*..********..**..****...**..*" : (66, 6.6), ".*.....*............*........." : (65, 6.5), ".*..*..*...*..*...*.**.......*" : (57, 5.7), ".*********************.*.*****" : (29, 2.9), "...*.............*.....*.....*" : (27, 2.7), ".*..*..*............*........*" : (24, 2.4), "....*......*......*..........." : (23, 2.3), ".*..**.*.*.*......*.*........*" : (19, 1.9), "................*..........*.." : (19, 1.9), ".*.....*...*......*.*........." : (14, 1.4), ".....*...*....*..............." : (10, 1.0), ".*..**.*.*.*..*...****...*...." : (8, 0.8), ".............*.............*.." : (8, 0.8), "...*...................*.....*" : (8, 0.8), ".*...*******..**..*.*.....*..." : (7, 0.7), ".*....***.**...*..*.......*..." : (7, 0.7), ".*****.*.*.*..*..*****.*.*..*." : (6, 0.6), "....*......*......*.*........." : (6, 0.6),
                "......*.*.*....*...*.....**..." : (6, 0.6), ".*****.*.*.*..*..*****.*.*..**" : (6, 0.6), ".....*...*....*.....*........." : (6, 0.6), "............**..*.....*.*..*.." : (6, 0.6), ".*.....**..*...*..*.......*..." : (5, 0.5), ".*.*.*******..**.*****.*.**..*" : (5, 0.5), ".*..*..*...*..*...*.*........." : (5, 0.5), ".....**.***...**....*.....*..." : (5, 0.5), "..*.*................*........" : (4, 0.4), ".....**.***...**..........*..." : (4, 0.4), ".........*....*.....*........." : (4, 0.4), ".........*....*..............." : (4, 0.4), ".*.*********..**.*****.*.**..*" : (4, 0.4), "......*.*.*..................." : (4, 0.4), ".......*..............*.*....." : (4, 0.4), "...................*.....**..." : (4, 0.4), "...................*.*...*...." : (4, 0.4), ".......*..............*.*....*" : (3, 0.3), "....*..............*.....*...." : (3, 0.3), ".*****.*.*.*..*..**.**.*....*." : (3, 0.3), "....*................*........" : (3, 0.3), ".*..*..*...*......*..*........" : (3, 0.3),
                "...*.**.***...*..*..*..*.....*" : (3, 0.3), "...*.**.*.*......*.....*.....*" : (3, 0.3), ".*.....*...*......**.*...**..." : (3, 0.3), ".....*...*..........*........." : (3, 0.3), ".***************.*************" : (3, 0.3), ".......*...........*..*.***..*" : (3, 0.3), ".................*...........*" : (3, 0.3), ".*...*******..**..*.**....*..." : (3, 0.3), ".............*......*......**." : (2, 0.2), ".....*...*...................*" : (2, 0.2), ".******.*******..*..**.*...**." : (2, 0.2), ".**************..*.***********" : (2, 0.2), "..**.............*...*.*....**" : (2, 0.2), "...*.*...*....*..*..*..*.....*" : (2, 0.2), "..*.*.................*.*....." : (2, 0.2), ".......*..............*.*.*..*" : (2, 0.2), ".........*..*................." : (2, 0.2), ".*.*********..**.*****.*.**.**" : (2, 0.2), ".***************.*.***********" : (2, 0.2), ".********.**...*.***.******..*" : (2, 0.2), ".*.*.*******..*..*****.*.**..*" : (2, 0.2), ".******.*.*......*...*.*......" : (2, 0.2),
                ".....**...*..................." : (2, 0.2), "...........*......*.*........*" : (2, 0.2), ".**.*..*...*......**.**.***.*." : (2, 0.2), ".***********..**.**********..*" : (2, 0.2), ".*....***.**...*..*.*.....*..." : (2, 0.2), ".....**..**.*................." : (2, 0.2), "..*..................*........" : (2, 0.2), ".............**.....*......**." : (2, 0.2), ".*..********..**..*.*.....*..." : (2, 0.2), "..**..*.*.*....*.*.*...*.**.*." : (2, 0.2), ".***********..*************.**" : (2, 0.2), ".......*..............*......." : (2, 0.2), ".********.*......*.*.******..*" : (2, 0.2), "...*.*...........*.....*......" : (2, 0.2), ".....*...*...........*........" : (2, 0.2), ".********.**.....*.*.******..*" : (2, 0.2), "...*.**.*.*......*.....*......" : (2, 0.2), "...*.**..**.*....*.....*......" : (2, 0.2), ".*..*..*...*......*.**........" : (2, 0.2), ".............*.............**." : (2, 0.2), "...*.**..**.*..........*......" : (2, 0.2), ".......*...*......*..........." : (2, 0.2),
                "...*..*.*.*......*.....*.....*" : (2, 0.2), ".*.....*...*......*.......*..." : (2, 0.2), "...*.**...*......*.....*......" : (2, 0.2), "...*.*...........*.....*.....*" : (2, 0.2), ".***********..*..**********.**" : (2, 0.2), ".***.******.*****.*.***.******" : (1, 0.1), ".*....*****..*.**...**..***.**" : (1, 0.1), "...*.*............*...*......." : (1, 0.1), ".***.******..*.**.*.***.******" : (1, 0.1), ".*....*****..*.**...**..***.*." : (1, 0.1), ".......*.....*............*..." : (1, 0.1), ".......*...*......*.........*." : (1, 0.1), "......*.*.*....*.............." : (1, 0.1), ".........*..........*........." : (1, 0.1), "....................*....*...." : (1, 0.1), ".*.........*.................." : (1, 0.1), ".......*.....*..*...*...***..." : (1, 0.1), ".......*.....*................" : (1, 0.1), "...........*.......*...*......" : (1, 0.1), "...........*.......*.........." : (1, 0.1), ".......*.....*......*....**..." : (1, 0.1), ".*.*********..**.**********..*" : (1, 0.1),
                ".**.*..............*.**.***..." : (1, 0.1), ".***.******..****.*.***.******" : (1, 0.1), "...*..............*...*......." : (1, 0.1), ".......*...*.................." : (1, 0.1), "....*.................*.*....." : (1, 0.1), "..*****.***.***..*..**.*...**." : (1, 0.1), "....*............*............" : (1, 0.1), ".*....***.**......*.......*..." : (1, 0.1), ".........*...........*........" : (1, 0.1), "...*.**..**.***..*..*..*...**." : (1, 0.1), ".....**...*.*................." : (1, 0.1), "..*..................*......*." : (1, 0.1), ".*.....*...*......**.....**..." : (1, 0.1), ".**********.*******.***.******" : (1, 0.1), "..***............*.*.*.*.*..**" : (1, 0.1), "..*...................*.*....." : (1, 0.1), ".*....*.*.*....*............*." : (1, 0.1), ".*.*.******..*.**.*.***.***.**" : (1, 0.1), ".*....*...*....*............*." : (1, 0.1), ".*....*...*....*.............." : (1, 0.1), ".....**..**.*..........*......" : (1, 0.1), ".*...*******...*..*.*.....*..." : (1, 0.1),
                "...................*.........*" : (1, 0.1), ".******.*****....*...*.*......" : (1, 0.1), ".......*.....*......*...***..." : (1, 0.1), "...*.**.***.*....*.....*......" : (1, 0.1), ".*.*.*.*.*.*..*..*****.*.**..*" : (1, 0.1), "......*.*.*...**..........*..." : (1, 0.1), "..*.*......*.........*........" : (1, 0.1), ".....*........*..............." : (1, 0.1), ".*..**.*.*.*......*.**........" : (1, 0.1), "..*........................*.." : (1, 0.1), "......*.*.*..................*" : (1, 0.1), ".*........*....*.............." : (1, 0.1), ".************..*.**********..*" : (1, 0.1), ".************.*************.**" : (1, 0.1), ".********.**...*.*.*.******..*" : (1, 0.1), ".*.................*.*...**..." : (1, 0.1), ".*.................*.....**..." : (1, 0.1), ".********.*......*...****.*..*" : (1, 0.1), ".*...*******..**..****...**..." : (1, 0.1), "..*****.*.*......*...*.*......" : (1, 0.1), ".*....***.**..**..*.......*..." : (1, 0.1), ".*....**..**......*..........." : (1, 0.1),
                "...*.**.***.***..*..*..*...**." : (1, 0.1), ".*..**.*.*.*..*...*.*........." : (1, 0.1), ".************..*.**********.**" : (1, 0.1), ".***.*******..**.**.**.*..*.**" : (1, 0.1), ".........*..*.......*........." : (1, 0.1), ".............**............*.." : (1, 0.1), ".***********.*****************" : (1, 0.1), ".***********.***.*************" : (1, 0.1), ".********.**.*.*.***.*********" : (1, 0.1), ".********.**...*.***.******.**" : (1, 0.1), ".********.**.....***.******..*" : (1, 0.1), ".*.*.**.*.*......*.....*......" : (1, 0.1), ".*....*.***....*.....*......*." : (1, 0.1), ".***********..**.*****.*.**..*" : (1, 0.1), ".......*...........*..*.*.*..*" : (1, 0.1), ".......*...........*..*.*....*" : (1, 0.1), ".......*...........*..*......*" : (1, 0.1), ".*.............*.............." : (1, 0.1), ".********.**.....**********.**" : (1, 0.1), ".********.**.....***.******.**" : (1, 0.1), ".*.*.**.***.*....*.....*......" : (1, 0.1), ".**.*..*...*......**.**.***..." : (1, 0.1),
                ".**.*..*...*.......*.**.***..." : (1, 0.1), ".**.*..*...........*.**.***..." : (1, 0.1), ".**.*..............*.*...**..." : (1, 0.1), ".**.*................*........" : (1, 0.1), ".*..*................*........" : (1, 0.1), "..................*...*......." : (1, 0.1), ".*.*.*******..**.**.*..*..*..*" : (1, 0.1), "....*...............*........*" : (1, 0.1), "........*....................*" : (1, 0.1), ".*.....**..*......*.......*..." : (1, 0.1),                     
                }),
            )
        for i in test_cases:
            self.check_group_freqs(i[0], i[1])            

if __name__ == "__main__":
    unittest.main()
    
    
                