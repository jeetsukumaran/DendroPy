#! /usr/bin/env python

############################################################################
##  parse_trees.py
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
Tree parser performance evaluation.
"""

from optparse import OptionParser
from optparse import OptionGroup
import sys
import os
import datetime

from dendropy import nexus


def main():   
    usage = "%prog [options] <TREES FILE> [<TREES FILE> [<TREES FILE> [...]]"
    parser = OptionParser(usage=usage, add_help_option=True)
    parser.add_option('-v', '--verbose', 
                      action='store_true', 
                      dest='verbose',
                      default=False,
                      help="suppress progress messages")     

    (opts, args) = parser.parse_args()
    
    support_filepaths = []        
    missing = False 
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))        
        if not os.path.exists(fpath):
            sys.stderr.write('File not found: "%s".\n' % fpath)
            missing = True
        else:
            support_filepaths.append(fpath)

    if len(support_filepaths) == 0:
        sys.stderr.write("No valid tree files specified or found.")
        sys.exit(1)            

        
#     tree_source = MultiFileTreeIterator(sources=support_filepaths,
#                                         core_iterator=nexus.iterate_over_trees, 
#                                         progress_func=progress_func)

    start_time = datetime.datetime.now()
    total_trees_read = 0
    for tree_file in support_filepaths:
        f = open(tree_file, "rU")
        for idx, tree in enumerate(nexus.iterate_over_trees(f)):
            pass
        total_trees_read += idx            
        
    end_time = datetime.datetime.now()
    parse_time = end_time-start_time
    
    hours, mins, secs = str(end_time-start_time).split(":")
    parse_seconds = float(hours * 3600) + float(mins * 60) + float(secs)
    
    if not not opts.verbose:
        sys.stderr.write("---\n\n")
    
    sys.stdout.write("%s\t%s\n" % ((", ".join(support_filepaths)), parse_seconds))
    
    if not not opts.verbose:
        sys.stderr.write("\nTotal trees parsed: %d\n" % tree_source.total_trees_read)
        sys.stderr.write("          Began at: %s.\n" % (start_time.isoformat(' ')))
        sys.stderr.write("          Ended at: %s.\n" % (end_time.isoformat(' ')))
        sys.stderr.write("        Parse time: %s hour(s), %s minute(s), %s second(s).\n" % (hours, mins, secs))
        sys.stderr.write("                   [= %s seconds]\n" % parse_seconds)
        
if __name__ == "__main__":
    main()

