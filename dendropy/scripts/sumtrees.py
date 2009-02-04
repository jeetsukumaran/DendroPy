#! /usr/bin/env python

############################################################################
##  sumtrees.py
##
##  Copyright 2008 Jeet Sukumaran.
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
CLI wrapper for tree summarization.
"""

__DEBUG__ = True

import os
import sys
import textwrap
from optparse import OptionParser
from optparse import OptionGroup

import datetime
import time
import socket
try:
    import getpass
except:
    pass
import platform    
    
import dendropy    
from dendropy import nexus
from dendropy import splits
from dendropy import treesum
from dendropy import datasets
from dendropy import trees

_program_name = 'SumTrees'
_program_subtitle = 'Phylogenetic Tree Split Support Summarization'
_program_date = 'Jan 12 2009'
_program_version = 'Version 1.0.5 (%s)' % _program_date
_program_author = 'Jeet Sukumaran'
_program_contact = 'jeetsukumaran@gmail.com'
_program_copyright = "Copyright (C) 2008 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

def show_splash(dest=sys.stderr, extended=False):
    lines = []
    lines.append("%s - %s" % (_program_name, _program_subtitle))
    lines.append("%s" % _program_version)
    lines.append("By %s" % _program_author)
    lines.append("(using the DendroPy Phylogenetic Computation Library Version %s)" % (dendropy.PACKAGE_VERSION))
    if extended:
        lines.append('')
        lines.extend(_program_copyright.split('\n'))
    header_max = max([len(i) for i in lines]) + 1
    sbars = '=' * header_max
    dest.write("%s\n" % sbars)
    dest.write("%s\n" % ('\n'.join(lines)))
    dest.write("%s\n\n" % sbars)       
        
class Messenger(object):

    def __init__(self, quiet=False, dest1=sys.stderr, dest2=None):
        self.quiet = quiet
        self.dest1 = dest1
        self.dest2 = dest2

    def send_multi(self, msg, wrap=0, newline=True, force=False):
        for line in msg:
            self.send(msg=line, wrap=wrap, newline=newline, force=force)
            
    def send(self, msg, wrap=0, newline=True, force=False):                
        if wrap:
            msg = textwrap.fill(msg, width=70)
        if newline:
            suffix = "\n"
        else:
            suffix = ""           
        if force or not self.quiet:            
            if self.dest1:
                self.dest1.write(msg + suffix)
        if self.dest2:
            self.dest2.write(msg + suffix)
            
    def send_formatted(self, msg, force=False):
        self.send(msg, wrap=True, force=force)
        
    def send_error(self, msg, wrap=False):
        self.send(msg, wrap=wrap, force=True)  
        
def main_cli():
    
    description =  '%s %s %s' % (_program_name, _program_version, _program_subtitle)    
    usage = "%prog [options] <TREES FILE> [<TREES FILE> [<TREES FILE> [...]]"
    
    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    sum_tree_optgroup = OptionGroup(parser, 'Summarization Options')    
    parser.add_option_group(sum_tree_optgroup)                      
    sum_tree_optgroup.add_option('-b', '--burnin', 
                        action='store',
                        dest='burnin',
                        type='int', # also 'float', 'string' etc.
                        default=0, 
                        help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]') 

    target_tree_optgroup = OptionGroup(parser, 'Target Tree Options')    
    parser.add_option_group(target_tree_optgroup)
    target_tree_optgroup.add_option('-t','--target',  
                  dest='target_tree_filepath',
                  default=None,
                  help="path to optional target, model or best topology tree file (Newick or NEXUS format) "
                       + "to which support will be mapped; " 
                       + "if not given, then a majority-rule clade consensus tree will be constructed based on the "
                       + "all the trees given in the support tree files (except for those discarded as burn-ins), "
                       + "and this will be used as the target tree")  
    target_tree_optgroup.add_option('-f', '--min-clade-freq', 
                      dest='min_clade_freq',
                      type='float', 
                      default=0.50,
                      metavar='#.##',
                      help="minimum frequency or probability for a clade or a split to be included in the consensus tree, if used [default=%default]") 
    target_tree_optgroup.add_option('--no-branch-lengths',  
                      action='store_true', 
                      dest='no_branch_lengths',
                      default=False,
                      help="by default, if using a consensus tree as the target tree, branch lengths will be the mean of the lengths " \
                          + "of the given branch across all trees considered; this option forces branch " \
                          + "lengths to be unspecified (obviously, this is only applicable if you do not ask the support to be mapped as "  \
                          + "branch lengths)")
    output_tree_optgroup = OptionGroup(parser, 'Output Tree Options')    
    parser.add_option_group(output_tree_optgroup)          
    output_tree_optgroup.add_option('-l','--support-as-labels',  
                      action='store_true', 
                      dest='support_as_labels',
                      default=True,
                      help="indicate branch support as internal node labels [default=%default]")            
    output_tree_optgroup.add_option('-v','--support-as-lengths',  
                      action='store_false', 
                      dest='support_as_labels',
                      default=True,
                      help="indicate branch support as branch lengths (otherwise support will be indicated by internal node labels)")   
    output_tree_optgroup.add_option('-p', '--percentages',  
                      action='store_true', 
                      dest='support_as_percentages',
                      default=False,
                      help="indicate branch support as percentages (otherwise, will report as proportions by default)")     
    output_tree_optgroup.add_option('-d', '--decimals', 
                      dest='support_label_decimals',
                      type='int', 
                      metavar='#',
                      default=2,
                      help="number of decimal places in indication of support values [default=%default]")  

                                            
    output_filepath_optgroup = OptionGroup(parser, 'Output File Options')    
    parser.add_option_group(output_filepath_optgroup)                      
    output_filepath_optgroup.add_option('-o','--output',  
                  dest='output_filepath',
                  default=None,
                  help="path to output file (if not given, will print to standard output)")                       
    output_filepath_optgroup.add_option('--no-taxa-block',  
                      action='store_false', 
                      dest='include_taxa_block',
                      default=True,
                      help="do not include a taxa block in the output treefile (otherwise will create taxa block by default)")      
    output_filepath_optgroup.add_option('--no-meta-comments',  
                      action='store_false', 
                      dest='include_meta_comments',
                      default=True,
                      help="include initial file comment annotating details of scoring operation")                      
    output_filepath_optgroup.add_option('-m', '--additional_comments',  
                      action='store', 
                      dest='additional_comments',
                      default=None,
                      help="additional comments to be added to the summary file")                                              
    output_filepath_optgroup.add_option('--newick', 
                      action='store_true', 
                      dest='phylip_format',
                      default=False,
                      help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")         
    output_filepath_optgroup.add_option('--phylip', 
                      action='store_true', 
                      dest='phylip_format',
                      default=False,
                      help="same as --newick")
    output_filepath_optgroup.add_option('-r', '--replace', 
                      action='store_true', 
                      dest='replace',
                      default=False,
                      help="replace/overwrite output file without asking if it already exists ")  
                                              
    run_optgroup = OptionGroup(parser, 'Program Run Options')    
    parser.add_option_group(run_optgroup)         
    run_optgroup.add_option('-q', '--quiet', 
                      action='store_true', 
                      dest='quiet',
                      default=False,
                      help="suppress progress messages") 
    run_optgroup.add_option('--ignore-missing-support', 
                      action='store_true', 
                      dest='ignore_missing_support',
                      default=False,
                      help="ignore missing support tree files (at least one must exist!)") 
    run_optgroup.add_option('--ignore-missing-target', 
                      action='store_true', 
                      dest='ignore_missing_target',
                      default=False,
                      help="ignore missing target tree file (will construct majority rule consensus tree if missing)") 
  
    (opts, args) = parser.parse_args()
    messenger = Messenger(quiet=opts.quiet)
    
    # splash 
    if not opts.quiet:
        show_splash(dest=sys.stderr, extended=False)
                                    
    ###################################################
    # Support file idiot checking
        
    support_filepaths = []        
    missing = False 
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))        
        if not os.path.exists(fpath):
            messenger.send_error('Support file not found: "%s"' % fpath)
            missing = True
        else:
            support_filepaths.append(fpath)
    if missing:
        messenger.send("")
        if opts.ignore_missing_support:
            pass
        else:
            messenger.send_formatted('Terminating due to missing support files. '
                   + 'Use the "--ignore-missing-support" option to continue even '
                   + 'if some files are missing.', force=True)
            sys.exit(1)
    if len(support_filepaths) == 0:
        messenger.send_formatted("No sources of support specified or could be found. "
        + "Please provide the path to at least one (valid and existing) file "
        + "containing non-parametric or MCMC tree samples "
        + "to summarize.", force=True)
        sys.exit(1)
        
    support_file_objs = [open(f, "r") for f in support_filepaths]

    ###################################################
    # Lots of other idiot-checking ...
    
    # target tree
    if opts.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(opts.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            messenger.send_error('Target tree file not found: "%s"\n' % target_tree_filepath)
            if opts.ignore_missing_target:                
                if not opts.quiet:
                    messenger.send('Will construct and use majority-rule consensus tree instead.\n')
                target_tree_filepath = None
            else:
                sys.exit(1)
    else:
        target_tree_filepath = None
                    
    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if os.path.exists(output_fpath):           
            if opts.replace:
                overwrite = 'y'
            else:
                messenger.send_error('Output path already exists: "%s"' % output_fpath)
                overwrite = raw_input("Overwrite (y/N)? ")
                messenger.send('')
            if not overwrite.lower().startswith("y"):
                sys.exit(1)
        output_dest = open(output_fpath, 'w')                
                
    ###################################################
    # Main work begins here: Count the splits
    
    start_time = datetime.datetime.now()
    
    comments = []
    tsum = treesum.TreeSummarizer()
    tsum.burnin = opts.burnin 
    tsum.support_as_labels = opts.support_as_labels 
    tsum.support_as_percentages = opts.support_as_percentages
    if not opts.support_as_percentages and opts.support_label_decimals < 2:
        messenger.send_error("(WARNING: proportions require that support will be reported to at least 2 decimal places)")
        opts.support_label_decimals = 2
    tsum.support_label_decimals = opts.support_label_decimals
    tsum.ignore_node_ages = True # until a more efficient implementation is developed
    if opts.quiet:
        tsum.verbose = False
        tsum.write_message = None
    else:
        tsum.verbose = True
        tsum.write_message = sys.stderr.write
        tsum.progress_message_prefix = ""
        tsum.progress_message_suffix = "\n"

    messenger.send("### COUNTING SPLITS ###\n")                
    split_distribution = tsum.count_splits_on_trees(tree_files=support_filepaths, 
                                                    tree_iterator=nexus.iterate_over_trees) 
        
    report = []
    report.append("%d trees read from %d files." % (tsum.total_trees_read, len(support_filepaths)))
    report.append("%d trees from each file requested to be ignored for burn-in." % (opts.burnin))
    report.append("%d trees ignored in total." % (tsum.total_trees_ignored))    
    report.append("%d trees considered in total for split support assessment." % (tsum.total_trees_counted))
    report.append("%d unique taxa across all trees." % len(split_distribution.taxa_block))
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = split_distribution.splits_considered()
    report.append("%d unique splits out of %d total splits counted." % (num_unique_splits, num_splits))
    report.append("%d unique non-trivial splits out of %d total non-trivial splits counted." % (num_nt_unique_splits, num_nt_splits))
        
    comments.extend(report)
    messenger.send("---")
    messenger.send_multi(report)
    messenger.send("")
    
    ###################################################
    #  Target tree and mapping
    
    if opts.support_as_percentages:
        support_units = "Percentage"
    else:        
        support_units = "Proportion (frequency or probability)"        
    if opts.support_as_labels:
        support_show = "node labels"
    else:
        support_show = "branch lengths"
    support_indication = "%s of support for each split indicated by %s" % (support_units, support_show)      
    
    tt_trees = []
    if target_tree_filepath is not None:
        messenger.send("### MAPPING SUPPORT TO TARGET TREE(S) ###\n")         
        tt_dataset = nexus.read_dataset(open(target_tree_filepath, 'r'))
        for tree_block in tt_dataset.trees_blocks:
            for tree in tree_block:
                tsum.map_split_support_to_tree(tree, split_distribution)
                tt_trees.append(tree)
        messenger.send('Parsed "%s": %d tree(s) in file' % (target_tree_filepath, len(tt_trees)))
        comments.append('Split support mapped to trees in:')
        comments.append('  - "%s" (%d trees)' % (os.path.abspath(target_tree_filepath), len(tt_trees)))
        comments.append(support_indication + ".")
    else:
        messenger.send("### CONSTRUCTING CLADE CONSENSUS TREE ###\n")
        if opts.min_clade_freq < 0.5:
            messenger.send("Minimum frequency threshold for clade inclusion is greater than 0.5: reset to 0.5.", force=True)
            min_freq = 0.5
        elif opts.min_clade_freq > 1.0:
            messenger.send("Maximum frequency threshold for clade inclusion is 1.0: reset to 1.0.", force=True)
            min_freq = 1.0
        else:            
            min_freq = opts.min_clade_freq
        tt_trees.append(tsum.tree_from_splits(split_distribution, 
                                              min_freq=min_freq, 
                                              include_edge_lengths=not opts.no_branch_lengths))
        report = []
        report.append('Consensus tree (%f clade frequency threshold) constructed from splits.' % min_freq)
        report.append(support_indication + ".")
        messenger.send_multi(report)
        comments.extend(report)
    messenger.send("")
                
    end_time = datetime.datetime.now()        
    
   
    ###################################################
    #  RESULTS    
            
    messenger.send("### RESULTS ###\n")
        
    final_run_report = []    
    final_run_report.append("Began at: %s." % (start_time.isoformat(' ')))
    final_run_report.append("Ended at: %s." % (end_time.isoformat(' ')))
    hours, mins, secs = str(end_time-start_time).split(":")
    run_time = "Run time: %s hour(s), %s minute(s), %s second(s)." % (hours, mins, secs)
    final_run_report.append(run_time)
                                
#     if not opts.output_filepath:
#         messenger.send('\n\n>>>>>>>>>>')
    
    output_dataset = datasets.Dataset()    
    taxa_block = output_dataset.add_taxa_block(taxa_block=split_distribution.taxa_block)
    trees_block = trees.TreesBlock()
    trees_block.taxa_block = taxa_block
    for tree in tt_trees:
        trees_block.append(tree)
    trees_block = output_dataset.add_trees_block(trees_block=trees_block)
        
    if opts.phylip_format:
        newick_writer = nexus.NewickWriter()
        newick_writer.write_dataset(output_dataset, output_dest)
    else:
        nexus_writer = nexus.NexusWriter()
        if opts.include_taxa_block:
            nexus_writer.simple = False
        else:
            nexus_writer.simple = True 
        if opts.include_meta_comments:
            nexus_writer.comment = []
            try:
                username = getpass.getuser()
            except:
                username = "a user"
            nexus_writer.comment.append("%s %s by %s." % (_program_name, _program_version, _program_author))
            python_version = sys.version.replace("\n", "").replace("[", "(").replace("]",")")            
            nexus_writer.comment.append("Running under Python %s on %s." % (python_version, sys.platform))               
            nexus_writer.comment.append("Executed on %s by %s@%s." % (platform.node(),  username, socket.gethostname()))         
            nexus_writer.comment.append("Basis of split support:")
            for support_file in support_filepaths:
                nexus_writer.comment.append('  - "%s"' % os.path.abspath(support_file))            
            nexus_writer.comment.extend(final_run_report)
            nexus_writer.comment.extend(comments)
        if opts.additional_comments:
            nexus_writer.comment.append("\n")
            nexus_writer.comment.append(opts.additional_comments)
            
        nexus_writer.write_dataset(output_dataset, output_dest)

    if not opts.output_filepath:
        #messenger.send('<<<<<<<<<')     
        pass
    else:
        messenger.send('Results written to: "%s".' % (output_fpath))
    messenger.send("")        
        
    ###################################################
    #  WRAP UP    
    messenger.send("### DONE ###\n")
    messenger.send_multi(final_run_report)        

if __name__ == '__main__':
    try:
        main_cli()
    except (KeyboardInterrupt, EOFError), e:
        sys.stderr.write("Terminating (user-abort).\n")
        sys.exit(1)
    except Exception, e:
        sys.stderr.write("Error encountered: %s : %s.\n" % (str(type(e)), e.message))
        raise # reraise exception, with correct traceback
