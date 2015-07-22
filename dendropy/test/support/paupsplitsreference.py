#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

import collections
import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
from dendropy.test.support import pathmap

_SPLITS_REFERENCE_FIELD_TYPES = [str, int, int, float, float]

def get_splits_reference(
        splits_filename,
        splits_dir=None,
        key_column_index=0):
    # Key columns are:
    #     0   : PAUP* bipartition string representation '....**...' etc.
    #     1   : unnormalized split bitmask (for rooted trees) == leafset_bitmask for all trees and split_bitmask for rooted trees
    #     2   : normalized split bitmask (for unrooted trees) == split_bitmask for unrooted trees
    #     3   : (weighted) counts
    #     4   : (weighted) frequencies
    if splits_dir is not None:
        splits_filepath = os.path.join(splits_dir, splits_filename)
    else:
        splits_filepath = pathmap.splits_source_path(splits_filename)
    d = collections.OrderedDict()
    with open(splits_filepath, "r") as src:
        for row in src:
            content = row.split("#")[0]
            if not content:
                continue
            fields = content.split("\t")
            assert len(fields) == 5, "{}: {}".format(content, fields)
            for idx, field in enumerate(fields):
                fields[idx] = _SPLITS_REFERENCE_FIELD_TYPES[idx](fields[idx])
            key = fields[key_column_index]
            d[key] = {
                "bipartition_string": fields[0],
                "unnormalized_split_bitmask": fields[1],
                "normalized_split_bitmask": fields[2],
                "count": fields[3],
                "frequency": fields[4]/100,
            }
    return d
