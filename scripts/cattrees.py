#! /usr/bin/env python

############################################################################
##  cattrees.py
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
Concatenate multiple tree files.
"""

import os
import sys
from optparse import OptionParser
from optparse import OptionGroup
import textwrap
    
import dendropy    
from dendropy import nexus
from dendropy import datasets
from dendropy import trees

_program_name = 'CatTrees'
_program_subtitle = 'Phylogenetic Tree File Concatenation'
_program_date = 'Oct 4 2008'
_program_version = 'Version 1.0.0 (%s)' % _program_date
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

    input_optgroup = OptionGroup(parser, 'Input File Options')    
    parser.add_option_group(input_optgroup)                         
    input_optgroup.add_option('-b', '--burnin', 
                        action='store',
                        dest='burnin',
                        type='int', # also 'float', 'string' etc.
                        default=0, 
                        help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]') 
                                            
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
                                            

                                                
    (opts, args) = parser.parse_args()
    messenger = Messenger(quiet=opts.quiet)
    
    # splash 
    if not opts.quiet:
        show_splash(dest=sys.stderr, extended=False)
                                    
    ###################################################
    # Tree file idiot checking
        
    tree_filepaths = []        
    missing = False 
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))        
        if not os.path.exists(fpath):
            messenger.send_error('Support file not found: "%s"' % fpath)
            missing = True
        else:
            tree_filepaths.append(fpath)
    if missing:
        messenger.send("")
        if opts.ignore_missing_support:
            pass
        else:
            messenger.send_formatted('Terminating due to missing support files. '
                   + 'Use the "--ignore-missing-support" option to continue even '
                   + 'if some files are missing.', force=True)
            sys.exit(1)
    if len(tree_filepaths) == 0:
        messenger.send_formatted("No sources of support specified or could be found. "
        + "Please specify path(s) to tree files to concatenate.", force=True)
        sys.exit(1)
        
    tree_file_objs = [open(f, "r") for f in tree_filepaths]

    ###################################################
    # Other prepping...
    
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
    # Main work begins here
    
    report = []
    trees_block = trees.TreesBlock()
    for tree_filepath_idx, tree_filepath in enumerate(tree_filepaths):
        messenger.send("Reading tree file %d of %d: %s" \
            % (tree_filepath_idx+1, len(tree_filepaths), tree_filepath))
        dataset = nexus.read_dataset(open(tree_filepath, "r"))
        trees_added = 0
        for tree_count, tree in enumerate(dataset.trees_blocks[0]):
            if tree_count >= opts.burnin:
                trees_block.append(tree)
                trees_added += 1
        message = "%d of %d trees added to tree collection (%d trees total)" \
            % (trees_added, len(dataset.trees_blocks[0]), len(trees_block))
        report.append("  - %s: %s" % (tree_filepath, message))
        messenger.send(" " + message)
        
    messenger.send("Collating trees ...")
    trees_block.normalize_taxa()
    output_dataset = datasets.Dataset()
    output_dataset.add_trees_block(trees_block=trees_block, taxa_block=trees_block.taxa_block)
    
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
            nexus_writer.comment.append("%s %s by %s." % (_program_name, _program_version, _program_author))
            nexus_writer.comment.append("Source trees:")
            nexus_writer.comment.extend(report)
        if opts.additional_comments:
            nexus_writer.comment.append("\n")
            nexus_writer.comment.append(opts.additional_comments)
            
        nexus_writer.write_dataset(output_dataset, output_dest)    

if __name__ == '__main__':
    try:
        main_cli()
    except (KeyboardInterrupt, EOFError), e:
        sys.stderr.write("Terminating (user-abort).\n")
        sys.exit(1)
    except Exception, e:
        sys.stderr.write("Error encountered: %s : %s.\n" % (str(type(e)), e.message))
        raise # reraise exception, with correct traceback
