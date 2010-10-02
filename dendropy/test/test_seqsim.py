#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Molecular character evolution tests.
"""

import unittest
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
from dendropy.test.support import runlevel
from dendropy.interop import paup
from dendropy import seqsim
import dendropy

if not paup.DENDROPY_PAUP_INTEROPERABILITY:
    _LOG.warn("PAUP interoperability not available: skipping sequence simulation tests")
else:

    tree_model_string = """
    #NEXUS
    BEGIN TAXA;
        DIMENSIONS NTAX=5;
        TAXLABELS
            A
            B
            C
            D
            E
      ;
    END;
    begin trees;
        tree true=(A:0.25,(B:0.25,(C:0.25,(D:0.25,E:0.25):0.25):0.25):0.25):0.25;
    end;
    """

    if runlevel.is_test_enabled(runlevel.SLOW, _LOG, __name__, "skipping all sequence generation frequency checking"):

        class SeqSimTest(unittest.TestCase):

            def setUp(self):
                self.tree_model = dendropy.Tree.get_from_string(tree_model_string, schema="NEXUS")
                assert self.tree_model.taxon_set is not None

            def estimate_params(self,
                seq_len=10000,
                kappa=1.0,
                base_freqs=[0.25, 0.25, 0.25, 0.25],
                unequal_base_freqs=True,
                gamma_rates=False,
                prop_invar=False):

                output_ds = seqsim.generate_hky_dataset(seq_len,
                    tree_model=self.tree_model,
                    kappa=kappa,
                    base_freqs=base_freqs)
                self.tree_model.reindex_taxa(output_ds.char_matrices[0].taxon_set)

                est_tree, mle = paup.estimate_model(char_matrix=output_ds.char_matrices[0],
                                                    tree_model=self.tree_model,
                                                    num_states=2,
                                                    unequal_base_freqs=unequal_base_freqs,
                                                    gamma_rates=gamma_rates,
                                                    prop_invar=prop_invar,
                                                    tree_est_criterion="likelihood",
                                                    tree_user_brlens=True,
                                                    paup_path='paup')

                return mle

            def testCharGen(self):
                kappas = [float(i)/10 for i in range(10, 200, 20)]
                base_freq_sets = [
                    [0.25, 0.25, 0.25, 0.25],
                    [0.4, 0.1, 0.4, 0.1]
                ]
                seq_len = 1000
                for kappa in kappas:
                    for bf in base_freq_sets:
                        mle = self.estimate_params(kappa=kappa, base_freqs=bf)
                        kappa_comp = "True = %f, Estimated = %f" % (kappa, mle['kappa'])
                        #kappa_threshold = 10000.00 / float(seq_len)
                        kappa_threshold = kappa * 0.2
                        _LOG.info("Kappa: %s [error threshold: %f]" % (kappa_comp, kappa_threshold))
                        self.assertTrue(abs(kappa-mle['kappa']) < kappa_threshold, \
                                    "Estimate over threshold (%f): %s" % (kappa_threshold, kappa_comp))

        if __name__ == "__main__":
            #@ Mark, if the slow test skip condition test goes here,
            # it only gets checked when the module is executed directly; which
            # means the condition does not get checked when invoked by setuptools
            unittest.main()
