#! /usr/bin/env python

############################################################################
##  paup.py
##
##  Part of the DendroPy phylogenetic computation library.
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

import re
import subprocess

def bipartitions(data_filepath,
                 tree_filepath,
                 min_clade_freq=0.5,
                 burnin=0,
                 paup_path="paup"):
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
    set warnreset=yes;
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
    lines = stdout.split('\n')
    tax_labels = []
    bipartitions = []
    bipartition_freqs = {}
    bipartition_counts = {}
    bipartition_pattern = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')
    taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
    for line in lines:
        bp_match = bipartition_pattern.match(line)
        if bp_match:
            bipartitions.append(bp_match)
            bipartition_counts[bp_match.group(1)] = int(bp_match.group(2))
            bipartition_freqs[bp_match.group(1)] = float(bp_match.group(3))
        else:
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                tax_labels.append(ti_match.group(2).strip())                
    return tax_labels, bipartitions, bipartition_counts, bipartition_freqs


bipartitions("dendropy/tests/data/anolis.chars.nexus", "dendropy/tests/data/anolis.mbcon.trees.nexus")
    