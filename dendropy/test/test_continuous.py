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

"""
Continuous character tests.
"""

import unittest
import inspect
import dendropy
from dendropy.test.support import pathmap
from dendropy.test.support.mockrandom import MockRandom
from dendropy.test.support import dendropytest
from dendropy.model import continuous

class BounceConstrainTest(unittest.TestCase):

    def runTest(self):
        output = [continuous._bounce_constrain(0.5, .5 - 0.05*i, 0.0, 1.0) for i in range(100)]
        ref = [(0.5, 0.5), (0.45, 0.47499999999999998), (0.4, 0.45), (0.34999999999999998, 0.42499999999999999), (0.29999999999999999, 0.4), (0.25, 0.375), (0.19999999999999996, 0.34999999999999998), (0.14999999999999997, 0.32499999999999996), (0.099999999999999978, 0.29999999999999999), (0.049999999999999989, 0.275), (0.0, 0.25), (0.05, 0.22954545454545455), (0.1, 0.21666666666666665), (0.15, 0.20961538461538459), (0.2, 0.20714285714285713), (0.25, 0.20833333333333331), (0.3, 0.2125), (0.35, 0.21911764705882353), (0.4, 0.2277777777777778), (0.45, 0.23815789473684212), (0.5, 0.25), (0.55, 0.2630952380952381), (0.6, 0.27727272727272728), (0.65, 0.29239130434782612), (0.7, 0.30833333333333335), (0.75, 0.32499999999999996), (0.8, 0.34230769230769231), (0.85, 0.36018518518518516), (0.9, 0.37857142857142867), (0.95, 0.39741379310344832), (1.0, 0.41666666666666669), (0.94999999999999996, 0.43467741935483872), (0.89999999999999991, 0.45), (0.84999999999999987, 0.46287878787878789), (0.79999999999999982, 0.47352941176470598), (0.75, 0.48214285714285721), (0.69999999999999996, 0.48888888888888893), (0.64999999999999991, 0.49391891891891893), (0.59999999999999987, 0.49736842105263157), (0.54999999999999982, 0.49935897435897436), (0.5, 0.5), (0.44999999999999973, 0.49939024390243897), (0.39999999999999991, 0.49761904761904763), (0.35, 0.49476744186046517), (0.29999999999999982, 0.49090909090909085), (0.25, 0.48611111111111116), (0.19999999999999973, 0.48043478260869565), (0.14999999999999991, 0.47393617021276602), (0.099999999999999645, 0.46666666666666667), (0.049999999999999822, 0.45867346938775505), (0.0, 0.45), (0.05, 0.4416666666666666), (0.1, 0.43461538461538463), (0.15, 0.42877358490566037), (0.2, 0.42407407407407405), (0.25, 0.42045454545454541), (0.3, 0.41785714285714287), (0.35, 0.41622807017543861), (0.4, 0.41551724137931034), (0.45, 0.41567796610169494), (0.5, 0.41666666666666669), (0.55, 0.41844262295081969), (0.6, 0.42096774193548392), (0.65, 0.4242063492063492), (0.7, 0.428125), (0.75, 0.43269230769230771), (0.8, 0.43787878787878792), (0.85, 0.44365671641791049), (0.9, 0.45), (0.95, 0.45688405797101456), (1.0, 0.4642857142857143), (0.94999999999999973, 0.47147887323943666), (0.89999999999999991, 0.4777777777777778), (0.84999999999999964, 0.48321917808219184), (0.79999999999999982, 0.48783783783783785), (0.75, 0.4916666666666667), (0.69999999999999973, 0.49473684210526325), (0.64999999999999991, 0.49707792207792206), (0.59999999999999964, 0.49871794871794872), (0.54999999999999982, 0.49968354430379752), (0.5, 0.5), (0.45, 0.49969135802469133), (0.39999999999999947, 0.49878048780487805), (0.34999999999999964, 0.49728915662650602), (0.29999999999999982, 0.49523809523809526), (0.25, 0.49264705882352944), (0.2, 0.48953488372093024), (0.14999999999999947, 0.48591954022988504), (0.099999999999999645, 0.4818181818181817), (0.049999999999999822, 0.47724719101123586), (0.0, 0.47222222222222221), (0.049999999999999822, 0.46730769230769231), (0.1, 0.4630434782608695), (0.15, 0.45940860215053758), (0.2, 0.45638297872340422), (0.25, 0.4539473684210526), (0.3, 0.45208333333333334), (0.35, 0.45077319587628867), (0.4, 0.45), (0.45, 0.44974747474747473)]
        assert_mat_approx_equal(output, ref, tester=self)

class KTBRatesTest(unittest.TestCase):

    def runTest(self):
        expected = [0.01, 0.00238612042332357, 0.00060046473606155989, 0.0032194690573412836, 0.0022181051458164411, 0.00017490090971284891, 0.022600788551531675, 0.014617722366202806, 0.00040617762499905189, 4.0911649525733303e-05, 0.51000000000000001, 0.12395537322536185, 0.41634777277597329, 0.018010492402013432, 0.077781889354748265, 0.039201420226780297, 0.0018067559608145154, 0.038064428026927236, 0.000468854154046173, 0.00027068933811426338, 1.01, 0.86242720407665974, 0.062893461668263165, 0.23334318118359176, 0.13647451390979612, 0.012908032460053776, 0.0036578960557443821, 0.024223063122821967, 0.14159784912709042, 0.082700071198721442, 1.51, 0.686445719853644, 3.3522062890772624, 2.4967837007875815, 0.056173355247824582, 0.0022005207631086566, 0.053495933645695409, 0.037267915554308549, 0.63252734941696775, 0.0078375963967802359, 2.0099999999999998, 2.1317309596278573, 1.5563076302280545, 0.83974983369284895, 0.023527204671999535, 0.022671617533422929, 0.13583576739840431, 0.33813396573401916, 0.036494567569726223, 0.00040501698299937316, 2.5099999999999998, 2.1818686975978907, 2.4488269166108831, 0.62237870951874186, 0.90954859339503202, 0.13834954150596895, 0.041786742411412607, 9.5704019633684995e-06, 0.19550604179155287, 2.3708404254904281, 3.0099999999999998, 2.4582221111663025, 0.32661845300155418, 0.059390570323355388, 0.9349562292466097, 0.45407795675810753, 0.27067855008521879, 0.1153677690757042, 0.085614112190463207, 0.029757425845535486, 3.5099999999999998, 5.5780589809075565, 0.14859621496133962, 1.1361638222601711, 0.47373053953720168, 12.559163966329541, 1.9442218373360938, 0.003988245875549439, 63.612522547967608, 0.62267918872625649, 4.0099999999999998, 10.110324482363417, 0.54728328439925211, 0.069082413746695395, 2.657686479590013, 0.16898912755533801, 0.39596036968888626, 0.030529950136951384, 0.0019733444343232381, 0.017929956714718788, 4.5099999999999998, 10.45969709101225, 3.3593154971792818, 1.0336291499040244, 0.84786240588642281, 2.278084793520756, 0.043356609471690816, 0.023168096061358393, 0.0059065469168308469, 0.62169733292424179]

        rng = MockRandom()
        y = [continuous._calc_KTB_rate(0.01+.5*i, 1, j, rng) for i in range(10) for j in range(10)]
        assert_vec_approx_equal(y, expected, tester=self)

        rng = MockRandom()
        y = [continuous._calc_KTB_rate(0.01+.5*i, j, 1, rng) for i in range(10) for j in range(10)]
        assert_vec_approx_equal(y, expected, tester=self)

        self.assertRaises(ValueError, continuous._calc_KTB_rate, 0, 1 , 1 , rng)

class TKPRatesTest(unittest.TestCase):

    def runTest(self):
        expected = [0.01, 0.0039340474963855641, 0.001632232380666595, 0.014428659286578023, 0.016389703355764434, 0.0021307292762983254, 0.45394897294494657, 0.48407248276349046, 0.022176546909804414, 0.0036827493270788596, 0.51000000000000001, 0.20436786045422725, 1.1317505850563239, 0.080717426949510368, 0.57473474392305157, 0.47757106516375003, 0.036289663562129378, 1.2605207376608101, 0.025598569446275875, 0.02436667769068299, 1.01, 1.4219020757916292, 0.17096215398172529, 1.0457715847484181, 1.008417839353775, 0.15725202748916964, 0.073470806288836357, 0.80215768313618252, 7.7309806110114048, 7.4444231676578108, 1.51, 1.1317576595037642, 9.1122414408448513, 11.189808222817929, 0.41506807319133626, 0.02680783090697382, 1.0744945504810126, 1.2341438671393068, 34.534823123534657, 0.70551794392946265, 2.0099999999999998, 3.5146301765484438, 4.2304827507410794, 3.7634976514794487, 0.17384383517242785, 0.27619684368030367, 2.7283343215701814, 11.197459097867025, 1.9925358755666318, 0.036458466937595801, 2.5099999999999998, 3.5972933316044275, 6.6566017084647564, 2.7893078600612546, 6.720705581299355, 1.6854424538625612, 0.83930915760416014, 0.00031692818644330357, 10.67426820212137, 213.41625387395686, 3.0099999999999998, 4.0529230827852576, 0.88784100563352941, 0.26617006989932601, 6.9084440279478647, 5.5318019658942186, 5.436724012051565, 3.8204558144079903, 4.6743721423293518, 2.678678109503108, 3.5099999999999998, 9.1966644910421707, 0.40392639090720356, 5.0919329843369292, 3.5004215324170675, 153.00193917129431, 39.050739500681885, 0.13207256469121023, 3473.1260500607032, 56.051794289673808, 4.0099999999999998, 16.669107027752826, 1.4876702070018706, 0.30960589864113691, 19.637794491060113, 2.0587090258674543, 7.9530766255051635, 1.0110130970614859, 0.1077409554922503, 1.6140032678015137, 4.5099999999999998, 17.245125079032153, 9.1315662720433046, 4.6324044639076911, 6.2649028812690863, 27.752754239037028, 0.87084078040787383, 0.76722197209428222, 0.32248653474293604, 55.963410447025701]

        rng = MockRandom()
        y = [continuous._calc_TKP_rate(0.01+.5*i, 1, j, rng) for i in range(10) for j in range(10)]
        assert_vec_approx_equal(y, expected, tester=self)

        rng = MockRandom()
        y = [continuous._calc_TKP_rate(0.01+.5*i, j, 1, rng) for i in range(10) for j in range(10)]
        assert_vec_approx_equal(y, expected, tester=self)

        try:
            continuous._calc_TKP_rate(0, 1, 1, rng)
        except ValueError:
            pass
        except OverflowError:
            pass

class KTBEvolveCrop(unittest.TestCase):

    def runTest(self):
        rng = MockRandom()
        newick = "((t5:1611.75,t6:1611.75):3922.93,((t4:1043.81,(t2:754.11,t1:754.11):2896.9):6584.0,t3:1702.21):3832.47);"
        tree = dendropy.Tree.get_from_string(newick, "newick")
        root = tree.seed_node
        root.mutation_rate = 1e-5
        root.mean_edge_rate = root.mutation_rate
        continuous.evolve_continuous_char(root, rng, roeotroe=0.01,
                            min_rate=1.0e-6, max_rate=1.0e-3, model='KTB',
                            time_attr='edge_length', val_attr='mutation_rate',
                            mean_val_attr='mean_edge_rate',
                            constrain_rate_mode="crop")
        for i in tree.preorder_node_iter():
            if i.edge_length is not None:
                i.edge_length *= i.mean_edge_rate

class KTBEvolveLinearBounce(unittest.TestCase):

    def runTest(self):
        rng = MockRandom()
        newick = "((t5:1611.75,t6:1611.75):3922.93,((t4:1043.81,(t2:754.11,t1:754.11):2896.9):6584.0,t3:1702.21):3832.47);"
        tree = dendropy.Tree.get_from_string(newick, "newick")
        root = tree.seed_node
        root.mutation_rate = 1e-5
        root.mean_edge_rate = root.mutation_rate
        continuous.evolve_continuous_char(root, rng, roeotroe=0.01,
                            min_rate=1.0e-6, max_rate=1.0e-3, model='KTB',
                            time_attr='edge_length', val_attr='mutation_rate',
                            mean_val_attr='mean_edge_rate',
                            constrain_rate_mode="linear_bounce")
        for i in tree.preorder_node_iter():
            if i.edge_length is not None:
                i.edge_length *= i.mean_edge_rate

class BifurcatingTreePICTest(dendropytest.ExtendedTestCase):

    def setUp(self):
        tree_str = "[&R] ((((Homo:0.21,Pongo:0.21)N1:0.28,Macaca:0.49)N2:0.13,Ateles:0.62)N3:0.38,Galago:1.00)N4:0.0;"
        data_str = """
    #NEXUS
    BEGIN DATA;
        DIMENSIONS  NTAX=5 NCHAR=2;
        FORMAT DATATYPE = CONTINUOUS GAP = - MISSING = ?;
        MATRIX
            Homo      4.09434   4.74493
            Pongo     3.61092   3.33220
            Macaca    2.37024   3.36730
            Ateles    2.02815   2.89037
            Galago   -1.46968   2.30259
        ;
    END;
    """
        taxa = dendropy.TaxonNamespace()
        self.tree = dendropy.Tree.get_from_string(tree_str, 'newick', taxon_namespace=taxa)
        self.char_matrix = dendropy.ContinuousCharacterMatrix.get_from_string(data_str,
                'nexus',
                taxon_namespace=taxa)
        self.pic = continuous.PhylogeneticIndependentConstrasts(tree=self.tree,
                char_matrix=self.char_matrix)
        self.expected_vals = []
        self.expected_vals.append({
            # state, corrected edge length, contrast, contrast_var
            "N1": (3.852630000, 0.385000000, 0.483420000, 0.420000000),
            "N2": (3.200378400, 0.345600000, 1.482390000, 0.875000000),
            "N3": (2.780823579, 0.601905551, 1.172228400, 0.965600000),
            "N4": (1.183724613, 0.375743470, 4.250503579, 1.601905551),
            })
        self.expected_vals.append({
            # state, corrected edge length, contrast, contrast_var
            "N1": (4.038565000, 0.385000000, 1.412730000, 0.420000000),
            "N2": (3.743208400, 0.345600000, 0.671265000, 0.875000000),
            "N3": (3.437967150, 0.601905551, 0.852838400, 0.965600000),
            "N4": (3.011356599, 0.375743470, 1.135377150, 1.601905551),
            })

    def testTreeValues(self):
        for cidx in range(self.char_matrix.vector_size):
            ctree = self.pic.contrasts_tree(character_index=cidx,
                    annotate_pic_statistics=True,
                    state_values_as_node_labels=False,
                    corrected_edge_lengths=False)
            for nd in ctree.postorder_internal_node_iter():
                vals = (nd.pic_state_value,
                        nd.pic_corrected_edge_length,
                        nd.pic_contrast_raw,
                        nd.pic_contrast_variance)
                exp_vals = self.expected_vals[cidx][nd.label]
                for vidx, val in enumerate(vals):
                    self.assertAlmostEqual(vals[vidx], exp_vals[vidx])

class MultifurcatingTreePICTest(dendropytest.ExtendedTestCase):

    def setUp(self):
        tree_str = "[&R] ((((Homo:0.21,Bogus1:0.23,Pongo:0.21)N1:0.28,Bogus2:0.49,Macaca:0.49)N2:0.13,Bogus3:0.62,Ateles:0.62)N3:0.38,Galago:1.00)N4:0.0;"
        data_str = """
    #NEXUS
    BEGIN DATA;
        DIMENSIONS  NTAX=8 NCHAR=2;
        FORMAT DATATYPE = CONTINUOUS GAP = - MISSING = ?;
        MATRIX
            Homo      4.09434   4.74493
            Pongo     3.61092   3.33220
            Macaca    2.37024   3.36730
            Ateles    2.02815   2.89037
            Galago   -1.46968   2.30259
            Bogus1    2.15      2.15
            Bogus2    2.15      2.15
            Bogus3    2.15      2.15
        ;
    END;
    """
        taxa = dendropy.TaxonNamespace()
        self.tree = dendropy.Tree.get_from_string(tree_str, 'newick', taxon_namespace=taxa)
        self.char_matrix = dendropy.ContinuousCharacterMatrix.get_from_string(data_str,
                'nexus',
                taxon_namespace=taxa)

    def testErrorOnDefault(self):
        pic = continuous.PhylogeneticIndependentConstrasts(tree=self.tree,
                char_matrix=self.char_matrix)
        self.assertRaises(ValueError, pic.contrasts_tree, 1)

    def testError(self):
        pic = continuous.PhylogeneticIndependentConstrasts(tree=self.tree,
                char_matrix=self.char_matrix,
                polytomy_strategy="error")
        self.assertRaises(ValueError, pic.contrasts_tree, 1)

    def testIgnore(self):
        pic = continuous.PhylogeneticIndependentConstrasts(tree=self.tree,
                char_matrix=self.char_matrix,
                polytomy_strategy="Ignore")
        ctree = pic.contrasts_tree(1)

    def testResolve(self):
        pic = continuous.PhylogeneticIndependentConstrasts(tree=self.tree,
                char_matrix=self.char_matrix,
                polytomy_strategy="Resolve")
        ctree = pic.contrasts_tree(1)

def approx_equal(x, y, tol=1e-5):
    "Returns True if x and y differ by less than tol"
    return (abs(x - y) < tol)

def vec_approx_equal(x, y, tol=1e-5):
    """Returns True if each element in the iterable ``x`` differs by less than
    ``tol`` from the corresponding element in ``y``
    """
    if len(x) != len(y):
        return False
    for i, j in zip(x, y):
        if abs(i - j) >= tol:
            return False
    return True

def mat_approx_equal(x, y, tol=1e-5):
    """Returns True if each cell in 2D iterable matrix ``x`` differs by less than
    ``tol`` from the corresponding element in ``y``
    """
    if len(x) != len(y):
        return False
    for row_x, row_y in zip(x, y):
        if len(row_x) != len(row_y):
            return False
        for i, j in zip(row_x, row_y):
            if abs(i - j) >= tol:
                return False
    return True

def _failure(tester, msg):
    calling_frame = inspect.currentframe().f_back.f_back
    co = calling_frame.f_code
    emsg = "%s\nCalled from file %s, line %d, in %s" % (msg, co.co_filename, calling_frame.f_lineno, co.co_name)
    if isinstance(tester, unittest.TestCase):
        tester.assertTrue(False, emsg)
    else:
        raise AssertionError(emsg)

def assert_approx_equal(x, y, tester=None, tol=1e-5):
    """Asserts that x and y are approximately equal.

    If ``tester`` is a unittest.TestCase then assertTrue is used; otherwise
    AssertionErrors are raised.
    """
    if abs(x - y) >= tol:
        _failure(tester, "%f != %f" % (x, y))

def assert_vec_approx_equal(x, y, tester=None, tol=1e-5):
    """Returns True if each element in the iterable ``x`` differs by less than
    ``tol`` from the corresponding element in ``y``

    If ``tester`` is a unittest.TestCase then assertTrue is used; otherwise
    AssertionErrors are raised.
    """
    if len(x) != len(y):
        _failure(tester, "vectors of numbers differ in length (%d vs %d)" % (len(x), len(y)))
    for n, itup in enumerate(zip(x, y)):
        i, j = itup
        if abs(i - j) >= tol:
            _failure(tester, "%f != %f at element %d" % (i, j, n))

def assert_mat_approx_equal(x, y, tester=None, tol=1e-5):
    """Returns True if each cell in 2D iterable matrix ``x`` differs by less than
    ``tol`` from the corresponding element in ``y``
    If ``tester`` is a unittest.TestCase then assertTrue is used; otherwise
    AssertionErrors are raised.
    """
    if len(x) != len(y):
        _failure(tester, "Matrices differs in length (%d vs %d)" % (len(x), len(y)))
    for n, row_tup in enumerate(zip(x, y)):
        row_x, row_y = row_tup
        if len(row_x) != len(row_y):
            _failure(tester, "row %d of matrix differs in length (%d vs %d)" % (n, len(row_x), len(row_y)))
        for col, cell_tup in enumerate(zip(row_x, row_y)):
            i, j = cell_tup
            if abs(i - j) >= tol:
                _failure(tester, "%f != %f for column %d of row %d" % (i, j, col, n))

if __name__ == "__main__":
    unittest.main()
