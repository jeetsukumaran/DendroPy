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
import csv

import unittest
import dendropy.tests
from dendropy import get_logger
_LOG = get_logger("PAUPWrapper")

from dendropy import datasets
from dendropy import dataio
from dendropy import taxa
from dendropy import splits
from dendropy import utils


if "PAUP_PATH" in os.environ:
    PAUP_PATH = os.environ["PAUP_PATH"]
else:
    PAUP_PATH = "paup"

    
###############################################################################
# HIGHER-LEVEL CONVENIENCE AND UTILITY METHODS 
    
def get_split_distribution(tree_filepaths, 
                           taxa_filepath, 
                           unrooted=True, 
                           burnin=0):
    """Returns a SplitDistribution object of splits calculated over
    specified trees"""
    p = Paup()
    p.stage_execute_file(taxa_filepath, clear_trees=True)      
    p.stage_list_taxa()        
    p.stage_load_trees(tree_filepaths=tree_filepaths, unrooted=unrooted, burnin=burnin)
    p.stage_count_splits()
    p.run()
    taxa_block = p.parse_taxa_block()
    tree_count, bipartition_counts = p.parse_group_freqs() 
    sd = build_split_distribution(bipartition_counts,
                                  tree_count,
                                  taxa_block,
                                  unrooted=unrooted)
    return sd                                  
    
###############################################################################
## PAUP* WRAPPER
        
class Paup(object):
    """ Wrapper around PAUP* """
    
    def __init__(self, paup_path=None):
        if paup_path is None:
            self.paup_path = PAUP_PATH
        else:
            self.paup_path = paup_path
        self.commands = []
        self.output = []
        
    ### WRAPPER OPERATIONS ###      
                    
    def run(self):
        """ executes list of commands in PAUP*, 
        return results of stdout """
        commands = "\n".join(self.commands) + "\n"
        paup_run = subprocess.Popen(['%s -n' % self.paup_path],
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        stdout, stderr = paup_run.communicate(commands)
        results = stdout.split("\n")
        if stderr:
            _LOG.error("\n*** ERROR FROM PAUP ***")
            _LOG.error(stderr)          
            _LOG.error("\n*** COMMANDS SENT TO PAUP ***\n")
            _LOG.error(commands)
            sys.exit(1)
        self.output.extend(results)
        return results
                       
    ### PAUP COMMANDS ###   
                 
    def stage_list_taxa(self):
        """ 
        Given a data file in memory, this gets PAUP* to print a list of 
        taxa that can be used to build a TaxaBlock later.
        """
        self.commands.append("[!TAXON LIST BEGIN]\ntstatus / full;\n[!TAXON LIST END]\n")
        
    def stage_count_splits(self, majrule_filepath=None, majrule_freq=0.5):
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
        paup_template.append("[!SPLITS COUNT BEGIN]")        
        paup_template.append("contree / strict=no %s showtree=no grpfreq=yes majrule=yes percent=%d;" \
            % (treefile, percent));
        paup_template.append("[!SPLITS COUNT END]")
        self.commands.extend(paup_template)
        
    def stage_execute_file(self,
                             filepath,
                             clear_trees=False):
        """Executes file, optionally clearing trees from file if requested"""
        self.commands.append("execute %s;" % filepath)
        if clear_trees:
            self.commands.append("cleartrees;")  
            
    def stage_deroot(self):
        self.commands.append("deroot;")          
    
    def stage_load_trees(self,
                           tree_filepaths,
                           unrooted=True,                           
                           burnin=0,
                           mode=7): # keep trees in memory, specify 3 to clear
        """
        Composes commands to load a set of trees into PAUP*, with the specified 
        number of burnin dropped. NOTE: Taxa Block must be active.
        """
        if isinstance(tree_filepaths, str):
            raise Exception("expecting list of filepaths, not string")
        if unrooted:
            rooting = "unrooted=yes"
        else:
            rooting = "rooted=yes"
        gettree_template = 'gett file= %%s storebrlens=yes warntree=no %s from=%d mode=%d;' % (rooting, burnin+1, mode)
        paup_template = []
        paup_template.append("set warnreset=no; set increase=auto; set warnroot=no;")
        for tree_filepath in tree_filepaths:
            paup_template.append(gettree_template % tree_filepath)                        
        self.commands.extend(paup_template)
        
    ### OUTPUT PARSERS ### 
                    
    def parse_taxa_block(self):
        """
        Given PAUP* output that includes a taxon listing as produced by
        `stage_list_taxa`, this parses out and returns a taxon block.
        """
        taxlabels = []
        taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
        idx = 0
        for line in self.output:
            idx += 1
            if line == "TAXON LIST BEGIN":
                break                  
        for line in self.output[idx:]:
            if line == "TAXON LIST END":
                break
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                taxlabels.append(ti_match.group(2).strip())                
        taxa_block = taxa.TaxaBlock() 
        for taxlabel in taxlabels:
            taxa_block.add_taxon(label=taxlabel)
        return taxa_block           
        
    def parse_group_freqs(self):
        """
        Given PAUP* output that includes a split counting procedure,
        this collects the splits and returns a dictionary of group strings and 
        their frequencies
        """        
        bipartitions = []
        bipartition_freqs = {}
        bipartition_counts = {}
        tree_count = None
        tree_count_pattern = re.compile('.*Majority-rule consensus of ([\d]*) tree.*', re.I)
        
        bipartition_section = re.compile('Bipartitions found in one or more trees and frequency of occurrence:')
        bp_full_row_with_perc_col = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')
        bp_full_row_with_no_perc_col = re.compile('([\.|\*]+)\s+([\d\.]+)')
        bp_row = re.compile('([\.|\*]+).*')

        # find tree count
        for idx, line in enumerate(self.output):
            tp_match = tree_count_pattern.match(line)
            if tp_match:
                break
        if not tp_match:
            raise Exception("Failed to find tree count in PAUP* output")
        tree_count = int(tp_match.group(1))
                
        while not bp_row.match(self.output[idx]):
            idx += 1
            
        split_idx = 0
        split_reps = {}
        for line in self.output[idx:]:
            if line == "SPLITS COUNT END":
                 break        
            bp_match = bp_full_row_with_perc_col.match(line)
            if not bp_match:
                bp_match = bp_full_row_with_no_perc_col.match(line)
            if bp_match:
                # full row, or end of partial rows      
                if len(split_reps) == 0:
                    split_rep = bp_match.group(1)
                else:
                    split_rep = split_reps[split_idx] + bp_match.group(1)
                bipartition_counts[split_rep] = int(bp_match.group(2))
                split_idx += 1
            else:
                # either (1) partial row or (2) break between sections
                bp_match = bp_row.match(line)
                if not bp_match:
                    split_idx = 0
                else:
                    if split_idx in split_reps:
                        split_reps[split_idx] += bp_match.group(1)
                    else:
                        split_reps[split_idx] = bp_match.group(1)
                    split_idx += 1                   
        return tree_count, bipartition_counts

###############################################################################
# UTILITY METHODS     
    
def build_split_distribution(bipartition_counts,
                             tree_count,
                             taxa_block,
                             unrooted=True):
    """
    Returns a populated SplitDistribution object based on the given 
    bipartition info.
    """
    sd = splits.SplitDistribution(taxa_block=taxa_block)
    sd.unrooted = unrooted
    sd.total_trees_counted = tree_count
    for g in bipartition_counts:
        sd.add_split_count(paup_group_to_mask(g, normalized=unrooted),
            bipartition_counts[g])
    return sd    
    
def paup_group_to_mask(group_string, normalized=False):
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
    clade_mask = int(group_string.replace("*", "1").replace(".", "0"), 2)
    if normalized:
        mask=((2 ** len(group_string)) -1)
        return utils.NormalizedBitmaskDict.normalize(clade_mask, mask)
    else:
        return clade_mask

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
    results = stdout.split('\n')
    tax_labels = []
    bipartitions = []
    bipartition_freqs = {}
    bipartition_counts = {}
    bipartition_pattern = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')
    bipartition_pattern2 = re.compile('([\.|\*]+)\s+([\d\.]+)')
    taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
    for line in results:
        bp_match = bipartition_pattern.match(line)
        if bp_match:
            bipartitions.append(bp_match.group(1))
            bipartition_counts[bp_match.group(1)] = int(bp_match.group(2))
            bipartition_freqs[bp_match.group(1)] = float(bp_match.group(3))
        else:
            bp_match2 = bipartition_pattern.match(line)
            if bp_match2:
                bipartitions.append(bp_match2.group(1))
                bipartition_counts[bp_match2.group(1)] = int(bp_match2.group(2))
                bipartition_freqs[bp_match2.group(1)] = float(bp_match2.group(2))                
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

class PaupWrapperDumbTests(unittest.TestCase):
    """
    Checks basic running of PAUP commands, and correct parsing/extraction of 
    output.
    """
    
    class SplitsTestCase(object):
    
        def __init__(self, tree_filepath, taxa_filepath, splitscsv_filepath, num_trees):		
            self.tree_filepath = tree_filepath
            self.taxa_filepath = taxa_filepath
            self.splitscsv_filepath = splitscsv_filepath
            self.num_trees = num_trees
            
        def _get_treefilepath(self):
            return self.__tree_filepath                  
        def _set_treefilepath(self, f):
            self.__tree_filepath = dendropy.tests.data_source_path(f)                 
        tree_filepath = property(_get_treefilepath, _set_treefilepath)            
            
        def _get_taxafilepath(self):
            return self.__taxa_filepath                  
        def _set_taxafilepath(self, f):
            if f is not None:
                self.__taxa_filepath = dendropy.tests.data_source_path(f)                 
            else:
                self.__taxa_filepath = None
        taxa_filepath = property(_get_taxafilepath, _set_taxafilepath)        
        
        def _get_splitscsvfilepath(self):
            return self.__splitscsv_filepath                  
        def _set_splitscsvfilepath(self, f):
            self.__splitscsv_filepath = dendropy.tests.data_source_path(f)                 
        splitscsv_filepath = property(_get_splitscsvfilepath, _set_splitscsvfilepath)        
    
    def setUp(self):
    
        self.splits_test_cases = [
            PaupWrapperDumbTests.SplitsTestCase("feb032009.tre", 
                "feb032009.tre",
                "feb032009.splits.csv", 
                100),   
            ]
        from dendropy.tests import is_test_enabled, TestLevel    
        if is_test_enabled(TestLevel.SLOW, _LOG, module_name=__name__, message="skipping large tree files"):
            self.splits_test_cases.append(
                PaupWrapperDumbTests.SplitsTestCase("terrarana.random.unrooted.100.tre", 
                    "terrarana.random.unrooted.100.tre", 
                    "terrarana.random.unrooted.100.splits.csv", 
                    100))
            self.splits_test_cases.append(
                PaupWrapperDumbTests.SplitsTestCase("anolis.mcmct.trees.nexus", 
                "anolis.chars.nexus", 
                "anolis.mcmct.trees.splits.csv", 
                1001)) 
           
        self.taxa_test_cases = (
            ("feb032009.tre",("T01", "T02", "T03", "T04", "T05", "T06",
            "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14",
            "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22",
            "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
            "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38",
            "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46",
            "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54",
            "T55", "T56", "T57", "T58", "T59")),
            ("primates.chars.nexus", ("Lemur catta", "Homo sapiens",
            "Pan", "Gorilla", "Pongo", "Hylobates", "Macaca fuscata",
            "Macaca mulatta", "Macaca fascicularis", "Macaca sylvanus",
            "Saimiri sciureus", "Tarsius syrichta", ))
        )               

    def check_taxa_block(self, filename, taxlabels):
        """Loads a taxa block from `filename`, make sure taxa returned match
        `taxlabels`"""
        p = Paup()
        p.stage_execute_file(dendropy.tests.data_source_path(filename))
        p.stage_list_taxa()
        p.run()
        taxa_block = p.parse_taxa_block()
        assert len(taxa_block) == len(taxlabels)
        for i, t in enumerate(taxa_block):
            assert t.label == taxlabels[i]
                        
    def check_group_freqs(self, treefile, taxafile, exp_tree_count, group_freqs, unrooted=True):
        """Calculates group frequencies from `filename`, make sure that 
        frequencies match `group_freqs` (given as dictionary of PAUP* group
        strings and their counts for the file)."""        
        p = Paup()
        p.stage_execute_file(taxafile, clear_trees=True)      
        p.stage_list_taxa()        
        p.stage_load_trees(tree_filepaths=[treefile], unrooted=unrooted)
        p.stage_count_splits()
        p.run()
        taxa_block = p.parse_taxa_block()
        tree_count, bipartition_counts = p.parse_group_freqs()

        assert exp_tree_count == tree_count, "%s != %s" % (exp_tree_count, tree_count)
        assert len(group_freqs) == len(bipartition_counts), "%d != %d" % (len(group_freqs), len(bipartition_counts))
        for g in group_freqs:
            assert g in bipartition_counts
            assert group_freqs[g] == bipartition_counts[g], \
                "%s != %s" % (group_freqs[g], bipartition_counts[g])
                
        sd = build_split_distribution(bipartition_counts,
                                      tree_count,
                                      taxa_block,
                                      unrooted=unrooted)
                                      
        sf = sd.split_frequencies                                      
        for g in bipartition_counts:
            s = paup_group_to_mask(g, normalized=unrooted)
            assert s in sd.splits
            assert s in sd.split_counts
            assert sd.split_counts[s] == bipartition_counts[g]
            assert sd.total_trees_counted == exp_tree_count
            self.assertAlmostEqual(sf[s], float(bipartition_counts[g]) / exp_tree_count)
        
    def testGroupRepToSplitMask(self):  
        for i in xrange(0xFF):       
            s = splits.split_as_string(i, 8, ".", "*")[::-1]
            r = paup_group_to_mask(s, normalized=False)
            assert r == i, "%s  =>  %s  =>  %s" \
                % (splits.split_as_string(i, 8), s, splits.split_as_string(r, 8))    
        for i in xrange(0xFF):     
            s = splits.split_as_string(i, 8, "*", ".")[::-1]
            r = paup_group_to_mask(s, normalized=True)
            normalized = utils.NormalizedBitmaskDict.normalize(i, 0xFF)
            assert r == normalized, "%s  =>  %s  =>  %s" \
                % (splits.split_as_string(i, 8), s, splits.split_as_string(normalized, 8))                 
        for i in xrange(0xFF):     
            s = splits.split_as_string(i, 8, ".", "*")[::-1]
            r = paup_group_to_mask(s, normalized=True)
            normalized = utils.NormalizedBitmaskDict.normalize(i, 0xFF)
            assert r == normalized, "%s  =>  %s  =>  %s" \
                % (splits.split_as_string(i, 8), s, splits.split_as_string(normalized, 8))            

    def testTaxaBlock(self):   
        for i in self.taxa_test_cases:
            self.check_taxa_block(i[0], i[1])

    def testGroupFreqs(self):
        for stc in self.splits_test_cases:
            splits = dict([ (s[0], int(s[1])) for s in csv.reader(open(stc.splitscsv_filepath, "rU"))])
            _LOG.info("Checking splits in %s" % stc.tree_filepath)
            self.check_group_freqs(stc.tree_filepath, 
                stc.taxa_filepath, stc.num_trees, splits)              

if __name__ == "__main__":
    unittest.main()
