#! /usr/bin/env python

############################################################################
##  compare-splits.py
##
##  Compares frequencies of splits across different tree files.
##
##  Copyright 2009 Jeet Sukumaran.
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
Constructs table of frequencies of splits in different empirical distributions
of trees.
"""

import sys
import os
from optparse import OptionGroup
from optparse import OptionParser

from dendropy.cli import confirm_overwrite, show_splash, Messenger
from dendropy import datasets
from dendropy import treesum
from dendropy import nexus

_prog_usage = '%prog [options] <tree-files>'
_prog_name = 'COMPARE-SPLITS'
_prog_subtitle = 'split frequencies comparison'
_prog_date = 'Jun 28 2009'
_prog_version = 'Version 1.0.0 (%s)' % _prog_date
_prog_description = 'constructs table of frequencies of splits in different empirical distributions of trees'
_prog_author = 'Jeet Sukumaran '
_prog_copyright = 'Copyright (C) 2009 Jeet Sukumaran.'

def main():
    """
    Main CLI handler.
    """
    
    parser = OptionParser(usage=_prog_usage, 
        add_help_option=True, 
        version=_prog_version, 
        description=_prog_description)
        
    parser.add_option('-b', '--burnin',
        action='store',
        dest='burnin',
        type=int,
        default=0,
        metavar='BURN-IN',
        help='number of initial trees to discard from each file (default=%default)')
        
    parser.add_option('-s','--separator',  
        dest='separator',
        default='\t',
        help="character to use to separate/delimit columns (default=<TAB>)")
        
    parser.add_option('-o','--output',  
        dest='output_filepath',
        default=None,
        help="path to output file (if not given, will print to standard output)")
        
    parser.add_option('--no-header-row',  
        dest='show_header_row',
        action="store_false",
        default=True,
        help="skip output of column headers")
        
    parser.add_option('--no-split-string',   
        dest='show_split_string',
        action="store_false",
        default=True,
        help="do not output split representation")    
        
    parser.add_option('--ignore-missing-files', 
        action='store_true', 
        dest='ignore_missing_files',
        default=False,
        help="ignore missing tree files (at least one must exist!)")
        
    parser.add_option('-r', '--replace', 
        action='store_true', 
        dest='replace',
        default=False,
        help="replace/overwrite output file without asking if it already exists ")
        
    parser.add_option('-q', '--quiet', 
        action='store_true', 
        dest='quiet',
        default=False,
        help="suppress progress messages")         
        
    (opts, args) = parser.parse_args()
    messenger = Messenger(quiet=opts.quiet)
    
    # splash 
    if not opts.quiet:
        show_splash(prog_name=_prog_name, 
            prog_subtitle=_prog_subtitle, 
            prog_version=_prog_version, 
            prog_author=_prog_author, 
            prog_copyright=_prog_copyright, 
            dest=sys.stderr, 
            extended=False)  
    
    missing = False 
    tree_filepaths = []
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))        
        if not os.path.exists(fpath):
            messenger.send_error(('Tree file not found: "%s"' % fpath))
            missing = True
        else:
            tree_filepaths.append(fpath)
    if missing:
        if opts.ignore_missing_files:
            pass
        else:
                messenger.send_formatted('Terminating due to missing tree files. '
                       + 'Use the "--ignore-missing-support" option to continue even '
                       + 'if some files are missing.', force=True)
                sys.exit(1)
        if len(support_filepaths) == 0:
            messenger.send_formatted("No sources of support specified or could be found. "
            + "Please provide the path to at least one (valid and existing) file "
            + "containing tree samples "
            + "to summarize.", force=True)
            sys.exit(1)
             
    tree_file_objs = [open(f, "rU") for f in tree_filepaths]  
    
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(output_fpath, messenger, opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)
        
    dataset = datasets.Dataset()
    taxa_block = dataset.add_taxa_block()
    split_dists = {}
    tsum = treesum.TreeSummarizer()
    splits_to_consider = set()
    for tfile_idx, tfile in enumerate(tree_file_objs):
        messenger.send("File %d of %d: %s" % (tfile_idx+1, len(tree_file_objs), tfile.name))
        tree_iterator = nexus.iterate_over_trees(tfile, taxa_block=taxa_block, from_index=opts.burnin+1)
        split_dists[tfile] = tsum.count_splits_on_trees(tree_iterator, split_distribution=None, trees_splits_encoded=False)
        split_dists[tfile].calc_freqs()
        splits_to_consider.update(split_dists[tfile].splits)          
        num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = split_dists[tfile].splits_considered()
        messenger.send("  Trees counted: %s" % split_dists[tfile].total_trees_counted)      
        messenger.send("  Total number of splits: %s" % num_splits)      
        messenger.send("  Total number of unique splits: %s" % num_unique_splits)      
        messenger.send("  Total number of non-trivial splits: %s" % num_nt_splits)
        messenger.send("  Total number of unique non-trivial splits: %s\n" % num_nt_unique_splits)        
        
    if opts.show_header_row:
        column_labels = [f.name for f in tree_file_objs]
        if opts.show_split_string:
            column_labels.insert(0, "Split")        
        output_dest.write(opts.separator.join(column_labels) + "\n")                
        
    for split in splits_to_consider:
        freqs = []
        for tfile in tree_file_objs:
            if split in split_dists[tfile].splits:
                freqs.append(split_dists[tfile].split_frequencies[split])
            else:
                freqs.append(0.0)
        row = [str(f) for f in freqs]
        if opts.show_split_string:
            row.insert(0, nexus.split_to_newick(split, taxa_block))
        output_dest.write(opts.separator.join(row) + "\n")           
                                    
if __name__ == "__main__":
    main()

    