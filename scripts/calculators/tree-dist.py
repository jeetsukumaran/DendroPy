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

from dendropy.utility.cli import confirm_overwrite
from dendropy.utility.messaging import ConsoleMessenger
from dendropy import DataSet
from dendropy.treecalc import symmetric_difference
import dendropy


_prog_usage = '%prog [options] <tree-files>'
_prog_name = 'tree-dist'
_prog_subtitle = 'command line calculation of tree-to-tree distances'
_prog_date = 'May 03, 2010'
_prog_version = 'Version 0.0.1 (%s)' % _prog_date
_prog_author = 'Jeet Sukumaran and Mark Holder'
_prog_description =  '%s %s %s' % (_prog_name, _prog_version, _prog_subtitle)
_prog_copyright = 'Copyright (C) 2009 Jeet Sukumaran and Mark Holder.'

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(usage=_prog_usage,
        add_help_option=True,
        version=_prog_version,
        description=_prog_description)
    parser.add_option('-f', '--file-format', dest='schema', default='newick',
        help='Either "nexus", "newick", or "nexml"')
    parser.add_option('-s','--separator',
        dest='separator',
        default='\t',
        help="character to use to separate/delimit columns (default=<TAB>)")

    parser.add_option('-o','--output',
        dest='output_filepath',
        default=None,
        help="path to output file (if not given, will print to standard output)")

    parser.add_option('--no-labels',
        dest='show_header_row',
        action="store_false",
        default=True,
        help="skip output of column and row labels")

    parser.add_option('-r', '--replace',
        action='store_true',
        dest='replace',
        default=False,
        help="replace/overwrite output file without asking if it already exists ")

    (opts, args) = parser.parse_args()
    messenger = ConsoleMessenger(quiet=False)
    schema = opts.schema.lower()

    tree_filepaths = []
    if not args:
        sys.exit('Expecting at least one file name as an argument')
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))
        if not os.path.exists(fpath):
            messenger.send_error(('Tree file not found: "%s"' % fpath))
            sys.exit(1)
        tree_filepaths.append(fpath)

    tree_file_objs = [open(f, "rU") for f in tree_filepaths]

    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if not confirm_overwrite(output_fpath, messenger, opts.replace):
            sys.exit(1)
        output_dest = open(output_fpath, "w")
            

    dataset = DataSet()
    ts = dendropy.TaxonSet()
    
    all_trees = []
    for tfile_idx, tfile in enumerate(tree_file_objs):
        curr_trees = dendropy.TreeList.get_from_stream(tfile, schema, taxon_set=ts)
        all_trees.extend(curr_trees)

    num_trees = len(all_trees)
    if opts.show_header_row:
        x = ['-'] + [str(i+1) for i in range(num_trees)]
        output_dest.write("%s\n" % opts.separator.join(x))


    td_mat = []
    for i in xrange(num_trees):
        if opts.show_header_row:
            output_dest.write(str(1+i) + opts.separator)
        td_row = []    
        for j in xrange(num_trees):
            if j == i:
               td_row.append('-')
            elif j < i:
                td_row.append(td_mat[j][i])
            else:
                d = symmetric_difference(all_trees[i], all_trees[j])
                td_row.append(str(d))
        td_mat.append(td_row)
        output_dest.write("%s\n" % opts.separator.join(td_row))
        

if __name__ == "__main__":
    main()


