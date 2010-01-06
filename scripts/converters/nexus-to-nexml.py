#! /usr/bin/env python

from dendropy import datasets
import os
import sys
from optparse import OptionGroup
from optparse import OptionParser

_prog_usage = '%prog <NEXUS-FILEPATH>'
_prog_version = 'NEXUS-to-NEXML Version 1.0'
_prog_description = 'reads a nexus/newick file and writes data in NeXML format to stdout'
_prog_author = 'Jeet Sukumaran'

def main():
    """
    Main CLI handler.
    """
    
    parser = OptionParser(usage=_prog_usage, 
        add_help_option=True, 
        version=_prog_version, 
        description=_prog_description)

    (opts, args) = parser.parse_args()
    
    if len(args) == 0:
        sys.stderr.write("Please specify a Newick/NEXUS file to convert.\n")
        sys.exit(1)
    fpath = os.path.expanduser(os.path.expandvars(args[0]))
    if not os.path.exists(fpath):
        sys.stderr.write('File not found: %s\n' % fpath)
        sys.exit(1)
    d = datasets.Dataset()
    d.read(open(fpath, "rU"), "nexus")
    d.write(sys.stdout, "nexml")
    
if __name__ == '__main__':
    main()

    
