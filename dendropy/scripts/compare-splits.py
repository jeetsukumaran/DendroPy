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
        
    parser.add_option('--no-headers',  
        dest='headers',
        action="store_false",
        default=True,
        help="skip output of column headers")
        
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
             
    tree_file_objs = [open(f, "r") for f in tree_filepaths]  
    
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(output_fpath, messenger, opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)    
    
if __name__ == "__main__":
    main()

    