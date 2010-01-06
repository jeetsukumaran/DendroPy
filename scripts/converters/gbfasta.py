#! /usr/bin/env python

############################################################################
##  gbfasta.py
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
Relabels sequences in a GenBank FASTA file.
"""

import re
import sys
from optparse import OptionGroup
from optparse import OptionParser

from dendropy import datasets

_prog_usage = '%prog [options] [<file>]'
_prog_version = 'GBFASTA Version 1.0'
_prog_description = 'Relabels sequences in GenBank FASTA file.'
_prog_author = 'Jeet Sukumaran'
_prog_copyright = 'Copyright (C) 2009 Jeet Sukumaran.'

def main():
    """
    Main CLI handler.
    """

    parser = OptionParser(usage=_prog_usage,
        add_help_option=True,
        version=_prog_version,
        description=_prog_description)

    parser.add_option('-n', '--nexus',
        action='store_const',
        dest='schema',
        const='NEXUS',
        default="NEXUS",
        help='output in NEXUS format (default)')

    parser.add_option('-p', '--phylip',
        action='store_const',
        dest='schema',
        const='PHYLIP',
        help='output in NEXUS format (default)')

    parser.add_option('-f', '--fasta',
        action='store_const',
        dest='schema',
        const='FASTA',
        help='output in FASTA format')

    (opts, args) = parser.parse_args()

    if len(args) == 0:
        sys.stderr.write("(reading from standard input)\n")
        input = sys.stdin
    else:
        input = open(args[0], "rU")

    output = sys.stdout

    fd = datasets.Dataset()
    fd.read(input, "DNAFASTA")
    pattern = re.compile("gi\|.+\|.+\|(.+)\|\S* ([\w\.]+) ([\w\.]+) (\w+).*")
    for t in fd.taxa_blocks[0]:
        m = pattern.match(t.label)
        t.label = m.groups(1)[1] + "_" + m.groups(1)[2] + "_" + m.groups(1)[0]

    fd.write(output, opts.schema)

if __name__ == '__main__':
    main()



