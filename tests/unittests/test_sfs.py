#! /usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
Tests of site freqruency spectrum calculation.
"""

import unittest
import os
import sys
import math
import dendropy
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap

class SinglePopSfsTests(unittest.TestCase):

    def read_expected_sfs(self, filename):
        filepath = pathmap.char_source_path(filename)
        with open(filepath) as src:
            return [int(v) for v in src.read().strip().split(",")]

    def test_single_pop_sfs(self):
        for test_data_name in (
                "sfs_test_single_pop_10x10",
                "sfs_test_single_pop_100x500_01",
                "sfs_test_single_pop_100x500_02",
                "sfs_test_single_pop_100x500_03",
                "sfs_test_single_pop_100x500_04",
                "sfs_test_single_pop_100x500_05",
                "sfs_test_single_pop_100x500_06",
                "sfs_test_single_pop_100x500_07",
                "sfs_test_single_pop_100x500_08",
                "sfs_test_single_pop_100x500_09",
                "sfs_test_single_pop_100x500_10",
                ):
            for data_type in ("dna", "std"):
                obs_data_path = pathmap.char_source_path(test_data_name + ".data.{}.fasta".format(data_type))
                if data_type == "dna":
                    obs_data = dendropy.DnaCharacterMatrix.get(path=obs_data_path, schema="fasta")
                else:
                    obs_data = dendropy.StandardCharacterMatrix.get(path=obs_data_path, schema="fasta")
                expected_folded_sfs = self.read_expected_sfs(test_data_name + ".sfs.folded.txt")
                obs_folded_sfs = obs_data.folded_site_frequency_spectrum(is_pad_vector_to_unfolded_length=True)
                self.assertEqual(obs_folded_sfs, expected_folded_sfs)
                k = int(math.ceil(len(obs_data)/2.0)) + 1
                expected_folded_sfs = expected_folded_sfs[:k]
                obs_folded_sfs = obs_data.folded_site_frequency_spectrum(is_pad_vector_to_unfolded_length=False)
                self.assertEqual(obs_folded_sfs, expected_folded_sfs)

if __name__ == "__main__":
    unittest.main()

