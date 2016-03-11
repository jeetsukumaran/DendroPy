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

import math
import unittest
import collections
import dendropy
from dendropy.calculate import treecompare
from dendropy.calculate.treecompare import TreeShapeKernel
from dendropy.calculate.treecompare import AssemblageInducedTreeManager
from dendropy.calculate.treecompare import AssemblageInducedTreeShapeKernel

class TreeShapeKernelBasicCalculationTest(unittest.TestCase):

    def test_small_single(self):
        ## Test based on:
        ##      KAMPHIR
        ##      By Art F.Y. Poon and Rosemary McCloskey
        ##      https://github.com/ArtPoon/kamphir.git
        T1_str = "[&R] ( ( A:0.5, B:0.25 )E:0.5, ( C:0.25, D:0.25 )F:0.5 )G;"
        T2_str = "[&R] ( ( ( A:0.25, B:0.25 )E:0.5, C:0.25 )F:0.5, D:0.25 )G;"
        taxon_namespace = dendropy.TaxonNamespace()
        T1 = dendropy.Tree.get_from_string(T1_str, "newick", taxon_namespace=taxon_namespace)
        T2 = dendropy.Tree.get_from_string(T2_str, "newick", taxon_namespace=taxon_namespace)
        tree_shape_kernel = TreeShapeKernel(decay_factor=0.5, gauss_factor=1)
        assert tree_shape_kernel(T1, T2) == 1.125 * (1+math.exp(-0.0625))

    def test_different_size_trees(self):
        # sigma = 1
        # gaussFactor = 1
        # decayFactor = 0.1
        trees_str = """\
        [&R] (T1:12.5150033896,((T2:9.56663932282,(((T3:1.07963880769,T4:1.07963880769):4.33934073323,(T5:0.649236454849,T6:0.649236454849):4.76974308607):1.9951772698,T7:7.41415681072):2.1524825121):0.294448382741,T8:9.86108770556):2.65391568409):1.88016636219;
        [&R] ((((T1:2.63300806671,(T2:1.36260395686,T3:1.36260395686):1.27040410984):13.0805656885,(((T4:6.49642565476,(T5:4.76391492099,(T6:2.31456992456,(T7:1.43347670305,(T8:0.845099591454,(T9:0.749944158281,T10:0.749944158281):0.0951554331727):0.588377111599):0.881093221511):2.44934499642):1.73251073377):1.78625805849,(T11:1.99279678602,T12:1.99279678602):6.28988692723):5.93122188014,T13:14.2139055934):1.49966816178):9.5804936484,(((((((T14:0.205658120519,T15:0.205658120519):7.75482186971,T16:7.96047999023):0.460360751119,(T17:2.05148090483,T18:2.05148090483):6.36935983652):7.37034613696,T19:15.7911868783):0.64452345349,T20:16.4357103318):3.12259410597,(T21:13.3918100022,T22:13.3918100022):6.16649443557):3.77060549196,(((((T23:0.117188712769,T24:0.117188712769):0.350994673936,T25:0.468183386705):0.434444802184,T26:0.902628188889):15.8023737799,T27:16.7050019688):5.11792644629,(T28:18.0216678715,T29:18.0216678715):3.80126054359):1.50598151462):1.96515747383):2.27498485797,T30:27.5690522615):2.7054878136;
        [&R] ((T1:18.4726052336,(T2:16.7820588123,T3:16.7820588123):1.69054642133):21.0361145203,(((T4:22.1454645535,T5:22.1454645535):1.47756990124,((((((T6:3.72910479068,T7:3.72910479068):0.635795000427,(T8:3.12137642206,T9:3.12137642206):1.24352336905):4.02980332076,T10:8.39470311187):0.544576038938,T11:8.93927915081):4.15398513837,((T12:1.29491896299,T13:1.29491896299):3.95016532419,(T14:3.55950898913,T15:3.55950898913):1.68557529805):7.848180002):6.79081967153,(T16:13.3961088469,((T17:7.5123548288,T18:7.5123548288):4.47081545633,(T19:5.42073528802,T20:5.42073528802):6.56243499711):1.41293856181):6.48797511377):3.73895049406):0.434178411336,(T21:6.05244405939,(T22:2.47879429586,T23:2.47879429586):3.57364976352):18.0047688067):15.4515068878):33.9171079393;
        [&R] ((((T1:12.2171017101,(T2:5.80709236615,T3:5.80709236615):6.4100093439):0.163771078375,(T4:1.48920874144,T5:1.48920874144):10.891664047):6.14111109308,((T6:18.0816707702,T7:18.0816707702):0.315080693117,(T8:9.85665216391,((T9:3.58056555459,T10:3.58056555459):0.704080460611,T11:4.2846460152):5.57200614872):8.54009929938):0.125232418213):0.445776242893,(((T12:10.3677145836,(T13:1.58010326492,T14:1.58010326492):8.78761131866):2.28511142688,T15:12.6528260105):1.47666702282,(T16:12.4604750725,(((T17:1.66562876613,T18:1.66562876613):2.60275247994,T19:4.26838124607):0.527073064007,(T20:0.771151073421,T21:0.771151073421):4.02430323665):7.66502076239):1.66901796081):4.83826709113):19.4421800373;
        [&R] (((T1:10.2887625813,((T2:3.40673675066,T3:3.40673675066):4.35878457692,T4:7.76552132758):2.52324125374):8.94208354363,((((((T5:1.25729509129,T6:1.25729509129):1.98815007493,T7:3.24544516622):9.92816453819,T8:13.1736097044):0.131783344535,((T9:5.12529436775,T10:5.12529436775):2.8254015761,T11:7.95069594385):5.3546971051):0.31354548999,((T12:6.23707496437,T13:6.23707496437):1.40621592147,T14:7.64329088584):5.97564765309):2.09487953589,(((T15:4.48342349888,T16:4.48342349888):1.33113258373,T17:5.81455608261):0.722465080174,(T18:3.26170665225,T19:3.26170665225):3.27531451053):9.17679691204):3.51702805012):12.2978992269,(T20:26.4312399558,((((((T21:0.115067790962,T22:0.115067790962):0.614386167869,T23:0.729453958831):10.0469766167,(((T24:5.54290468454,T25:5.54290468454):0.24279629634,((T26:0.0173560007622,T27:0.0173560007622):1.03837416324,T28:1.055730164):4.72997081688):0.935147415819,T29:6.7208483967):4.05558217886):5.34149070025,T30:16.1179212758):0.249537503423,((T31:7.68676243321,T32:7.68676243321):0.224359722764,T33:7.91112215597):8.45633662326):4.68102050585,(((T34:6.33634296608,T35:6.33634296608):11.7427817803,T36:18.0791247464):1.88140825859,(T37:3.77606556293,T38:3.77606556293):16.1844674421):1.08794628008):5.38276067077):5.097505396):12.8756560997;
        [&R] (((T1:21.6283980853,(T2:10.2523638386,T3:10.2523638386):11.3760342467):5.43593515808,(((T4:14.5815467291,((T5:10.2802760256,(T6:2.13975995776,T7:2.13975995776):8.14051606783):1.62385561861,T8:11.9041316442):2.67741508492):5.60705586175,(T9:5.28631833362,(T10:5.19018203656,((T11:0.654535739431,(T12:0.451288383503,T13:0.451288383503):0.203247355928):0.3240577623,T14:0.978593501732):4.21158853482):0.0961362970605):14.9022842572):6.30642464811,(T15:6.96742024513,T16:6.96742024513):19.5276069938):0.56930600442):9.39112113938,(((T17:9.86297817832,(T18:1.09548421955,(T19:0.163658446308,T20:0.163658446308):0.931825773239):8.76749395877):13.4127381553,(((((T21:0.162641367112,T22:0.162641367112):7.60612335657,(T23:1.19220280515,T24:1.19220280515):6.57656191854):4.64147681058,T25:12.4102415343):2.18305782204,(T26:13.4944558875,((T27:4.0278122205,(T28:3.83892525427,T29:3.83892525427):0.188886966234):6.19053229632,(T30:0.0720558467085,T31:0.0720558467085):10.1462886701):3.27611137066):1.09884346882):4.83967942511,(T32:13.108711714,(T33:1.21934185331,T34:1.21934185331):11.8893698607):6.32426706741):3.84273755221):12.2937144592,((T35:3.69731694953,(T36:1.08556862997,T37:1.08556862997):2.61174831956):2.01430134794,T38:5.71161829747):29.8578124954):0.886023589904):13.2943061023;
        [&R] (((((T1:2.69139513211,T2:2.69139513211):8.551485817,T3:11.2428809491):6.39974203014,((T4:17.0587824705,(((((T5:1.0038538602,(T6:0.65693286312,(T7:0.372810170449,T8:0.372810170449):0.284122692671):0.346920997081):0.992649847145,T9:1.99650370735):8.59545004131,T10:10.5919537487):3.93618536755,(T11:11.4771452793,(T12:11.3117342188,(((T13:1.44496178243,T14:1.44496178243):1.39593941106,T15:2.84090119349):7.36176225558,T16:10.2026634491):1.10907076976):0.16541106051):3.05099383686):0.446923840522,((T17:2.07495063443,(T18:1.03526596889,T19:1.03526596889):1.03968466554):3.02573649405,(T20:3.09572655531,T21:3.09572655531):2.00496057316):9.87437582825):2.08371951376):0.535710348712,T22:17.5944928192):0.0481301600476):0.485204130118,(T23:2.10622863693,T24:2.10622863693):16.0215984724):8.64518545213,((T25:12.6080734092,((T26:8.70251718564,T27:8.70251718564):1.05101612319,((T28:6.04918606932,((T29:3.85573365482,T30:3.85573365482):0.0757038480577,T31:3.93143750288):2.11774856645):2.49018043616,(T32:4.24643709587,T33:4.24643709587):4.29292940962):1.21416680334):2.85454010035):3.79696048295,((((T34:2.04492296218,T35:2.04492296218):1.53035052582,T36:3.57527348801):0.0643892487425,T37:3.63966273675):12.6279014452,T38:16.267564182):0.137469710183):10.3679786694):5.11271917475;
        [&R] ((((T1:7.07157972517,(T2:0.0218808463127,T3:0.0218808463127):7.04969887886):0.8655582016,(T4:3.02712438178,T5:3.02712438178):4.91001354499):1.21415733877,(((T6:0.0583489618229,T7:0.0583489618229):1.01650984984,T8:1.07485881166):4.07712736061,T9:5.15198617228):3.99930909326):7.26950003749,(((T10:0.424456263841,T11:0.424456263841):8.04936250082,T12:8.47381876466):3.5117568466,((T13:6.98368939693,(T14:0.777662333494,T15:0.777662333494):6.20602706344):2.06669325023,T16:9.05038264716):2.9351929641):4.43521969176):10.6934883223;
        [&R] (((T1:14.4454814308,((T2:3.56688601323,T3:3.56688601323):5.84751651584,T4:9.41440252907):5.03107890172):20.1795615891,(((((T5:4.15054680055,T6:4.15054680055):11.295603376,(T7:6.43442739414,T8:6.43442739414):9.01172278238):0.10084264253,(T9:8.4004080966,T10:8.4004080966):7.14658472244):2.73330732663,(((T11:6.10577944378,(T12:4.04456805007,(T13:2.80049127557,T14:2.80049127557):1.2440767745):2.06121139371):1.90192017764,(T15:7.79084487337,T16:7.79084487337):0.216854748047):7.40306903527,T17:15.4107686567):2.86953148899):11.2997061766,(T18:1.00480866062,T19:1.00480866062):28.5751976616):5.04503669769):0.0473375886983,((((T20:5.77874649173,(T21:5.42355087634,(T22:3.52279709572,(T23:2.90310692644,(T24:2.09879096434,T25:2.09879096434):0.804315962107):0.619690169282):1.90075378062):0.355195615385):4.90915591989,(T26:4.01715333576,T27:4.01715333576):6.67074907586):3.94107350362,(T28:12.4346731779,(T29:3.0690908858,T30:3.0690908858):9.36558229213):2.1943027373):12.9196714424,(T31:14.8264557361,(T32:0.0945338917738,T33:0.0945338917738):14.7319218443):12.7221916216):7.123733251):12.8703596545;
        """
        expected = [[ [] for i in range(10) ] for i in range(10)]
        expected[0][0] = 1.02943649816
        expected[0][1] = 0.629696880372
        expected[0][2] = 0.226436166785
        expected[0][3] = 0.504437998991
        expected[0][4] = 0.532944026681
        expected[0][5] = 1.06675009986
        expected[0][6] = 0.63997331036
        expected[0][7] = 0.533570048403
        expected[0][8] = 0.316006296153
        expected[0][9] = 0.763121914566
        expected[1][0] = 0.629696880372
        expected[1][1] = 5.73381840453
        expected[1][2] = 0.609230566852
        expected[1][3] = 1.36716715734
        expected[1][4] = 1.66142414246
        expected[1][5] = 3.3424585382
        expected[1][6] = 2.81647328252
        expected[1][7] = 1.6426813156
        expected[1][8] = 1.71534697686
        expected[1][9] = 2.11186097095
        expected[2][0] = 0.226436166785
        expected[2][1] = 0.609230566852
        expected[2][2] = 3.41707388515
        expected[2][3] = 0.965818487686
        expected[2][4] = 1.8148659483
        expected[2][5] = 1.09814847947
        expected[2][6] = 1.70646748651
        expected[2][7] = 0.543992787817
        expected[2][8] = 1.89068146553
        expected[2][9] = 1.2301930343
        expected[3][0] = 0.504437998991
        expected[3][1] = 1.36716715734
        expected[3][2] = 0.965818487686
        expected[3][3] = 3.52734267013
        expected[3][4] = 1.6121784729
        expected[3][5] = 1.9713400052
        expected[3][6] = 1.968844455
        expected[3][7] = 0.727319654436
        expected[3][8] = 1.07080235987
        expected[3][9] = 2.06317190944
        expected[4][0] = 0.532944026681
        expected[4][1] = 1.66142414246
        expected[4][2] = 1.8148659483
        expected[4][3] = 1.6121784729
        expected[4][4] = 7.13717804604
        expected[4][5] = 2.47204082769
        expected[4][6] = 2.17747252512
        expected[4][7] = 1.58536532878
        expected[4][8] = 2.31434377631
        expected[4][9] = 2.23786059877
        expected[5][0] = 1.06675009986
        expected[5][1] = 3.3424585382
        expected[5][2] = 1.09814847947
        expected[5][3] = 1.9713400052
        expected[5][4] = 2.47204082769
        expected[5][5] = 7.77913926686
        expected[5][6] = 2.87296840413
        expected[5][7] = 2.30072202059
        expected[5][8] = 2.12831286543
        expected[5][9] = 3.00420928696
        expected[6][0] = 0.63997331036
        expected[6][1] = 2.81647328252
        expected[6][2] = 1.70646748651
        expected[6][3] = 1.968844455
        expected[6][4] = 2.17747252512
        expected[6][5] = 2.87296840413
        expected[6][6] = 7.32395549398
        expected[6][7] = 1.3099686298
        expected[6][8] = 2.72462663407
        expected[6][9] = 2.60804726918
        expected[7][0] = 0.533570048403
        expected[7][1] = 1.6426813156
        expected[7][2] = 0.543992787817
        expected[7][3] = 0.727319654436
        expected[7][4] = 1.58536532878
        expected[7][5] = 2.30072202059
        expected[7][6] = 1.3099686298
        expected[7][7] = 3.00348381586
        expected[7][8] = 1.05294859724
        expected[7][9] = 1.30614199254
        expected[8][0] = 0.316006296153
        expected[8][1] = 1.71534697686
        expected[8][2] = 1.89068146553
        expected[8][3] = 1.07080235987
        expected[8][4] = 2.31434377631
        expected[8][5] = 2.12831286543
        expected[8][6] = 2.72462663407
        expected[8][7] = 1.05294859724
        expected[8][8] = 6.04463164133
        expected[8][9] = 2.07119628549
        expected[9][0] = 0.763121914566
        expected[9][1] = 2.11186097095
        expected[9][2] = 1.2301930343
        expected[9][3] = 2.06317190944
        expected[9][4] = 2.23786059877
        expected[9][5] = 3.00420928696
        expected[9][6] = 2.60804726918
        expected[9][7] = 1.30614199254
        expected[9][8] = 2.07119628549
        expected[9][9] = 4.7667655343
        trees = dendropy.TreeList.get(data=trees_str, schema="newick")
        tree_shape_kernel = TreeShapeKernel(
                sigma=1,
                gauss_factor=1,
                decay_factor=0.1,
                )
        for idx1, t1 in enumerate(trees):
            for idx2, t2 in enumerate(trees):
                self.assertAlmostEqual(tree_shape_kernel(t1, t2), expected[idx1][idx2])
                # print("{}, {} = {}".format(idx1+1, idx2+1, tree_shape_kernel(t1, t2)))

class AssemblageInducedTreeManagerTestBase(unittest.TestCase):

    GROUP_IDS = ("a", "b", "c", "d", "e")

    @staticmethod
    def get_random_tree():
        from dendropy.model import birthdeath
        tns = dendropy.TaxonNamespace()
        for group_id in AssemblageInducedTreeManagerTests.GROUP_IDS:
            for group_member in range(10):
                t = tns.require_taxon(label="{}{}".format(
                    group_id,
                    group_member))
        tree = dendropy.simulate.birth_death_tree(
                birth_rate=0.1,
                death_rate=0.0,
                taxon_namespace=tns)
        tree.assemblage_leaf_sets = []
        tree.assemblage_classification_regime_subtrees = []
        for group_id in AssemblageInducedTreeManagerTests.GROUP_IDS:
            node_filter_fn=lambda nd: nd.taxon is None or nd.taxon.label.startswith(group_id)
            subtree1 = tree.extract_tree(
                    node_filter_fn=node_filter_fn)
            assemblage_leaf_set = set()
            for leaf_nd in tree.leaf_node_iter():
                if node_filter_fn(leaf_nd):
                    assert leaf_nd.taxon.label.startswith(group_id)
                    assemblage_leaf_set.add(leaf_nd)
            tree.assemblage_leaf_sets.append(assemblage_leaf_set)
            tree.assemblage_classification_regime_subtrees.append(subtree1)
        assert len(tree.assemblage_classification_regime_subtrees) == len(AssemblageInducedTreeManagerTests.GROUP_IDS)
        return tree

    def validate_managed_trees(self, test_target, trees):
        self.assertEqual(test_target._num_assemblage_classifications, len(AssemblageInducedTreeManagerTests.GROUP_IDS))
        self.assertEqual(len(test_target._tree_assemblage_induced_trees_map), len(trees))
        for tree in trees:
            self.assertIn(tree, test_target._tree_assemblage_induced_trees_map)
            self.assertEqual(len(test_target._tree_assemblage_induced_trees_map[tree]), len(AssemblageInducedTreeManagerTests.GROUP_IDS))
            self.assertEqual(len(test_target._tree_assemblage_induced_trees_map[tree]), len(tree.assemblage_leaf_sets))
            induced_trees = test_target._tree_assemblage_induced_trees_map[tree]
            for ( induced_tree, group_id, original_leafset_nodes) in zip( induced_trees, AssemblageInducedTreeManagerTests.GROUP_IDS, tree.assemblage_leaf_sets):
                original_leafset = set(original_leafset_nodes)
                for leaf_nd in induced_tree.leaf_node_iter():
                    self.assertTrue(leaf_nd.taxon.label.startswith(group_id), leaf_nd.taxon.label)
                    original_node = leaf_nd.extraction_source
                    self.assertIn(original_node, original_leafset)
                    original_leafset.remove(original_node)
                self.assertEqual(len(original_leafset), 0)
                labels=[x.taxon.label for x in original_leafset_nodes]
                t2 = tree.extract_tree_with_taxa_labels(labels=labels)
                self.assertEqual(treecompare.weighted_robinson_foulds_distance(t2, induced_tree), 0.0)
                t3 = dendropy.Tree(tree)
                t3.retain_taxa_with_labels(labels=labels)
                # print(t3.as_string("newick"))
                # print(induced_tree.as_string("newick"))
                self.assertAlmostEqual(treecompare.weighted_robinson_foulds_distance(t3, induced_tree), 0.0)

class AssemblageInducedTreeManagerTests(AssemblageInducedTreeManagerTestBase):

    def test_basic_subtree_generation_and_caching(self):
        test_target = AssemblageInducedTreeManager()
        trees = []
        for idx in range(10):
            tree = self.get_random_tree()
            tree.induced_trees = test_target.generate_induced_trees(
                tree=tree,
                assemblage_leaf_sets=tree.assemblage_leaf_sets)
            self.assertEqual(len(tree.induced_trees), len(AssemblageInducedTreeManagerTests.GROUP_IDS))
            self.assertEqual(len(tree.induced_trees), len(tree.assemblage_leaf_sets))
            for induced_tree in tree.induced_trees:
                self.assertIn(induced_tree, test_target._tree_assemblage_induced_trees_map[tree])
            trees.append(tree)
        self.validate_managed_trees(
                test_target=test_target,
                trees=trees)

class AssemblageInducedTreeShapeKernelTests(AssemblageInducedTreeManagerTestBase):

    def validate_managed_and_cached_trees(self, test_target, trees):
        self.validate_managed_trees(
                test_target=test_target,
                trees=trees)
        for tree in trees:
            self.assertIn(tree, test_target._tree_cache)
            self.assertEqual(len(test_target._tree_assemblage_induced_trees_map[tree]), len(AssemblageInducedTreeManagerTests.GROUP_IDS))
            # for induced_tree in test_target._tree_assemblage_induced_trees_map[tree]:
            #     self.assertIn(induced_tree, test_target._tree_cache)

    def test_comparison_add_trees_directly(self):
        test_target = AssemblageInducedTreeShapeKernel()
        trees = [self.get_random_tree() for idx in range(10)]
        for tree in trees:
            test_target.update_assemblage_induced_tree_cache(
                    tree=tree,
                    assemblage_leaf_sets=tree.assemblage_leaf_sets)
            ## to check for inadvertent double adding
            test_target.update_assemblage_induced_tree_cache(
                    tree=tree,
                    assemblage_leaf_sets=tree.assemblage_leaf_sets)
        self.validate_managed_and_cached_trees(
                test_target=test_target,
                trees=trees)

    def test_comparison_add_trees_by_comparison(self):
        test_target = AssemblageInducedTreeShapeKernel()
        trees = [self.get_random_tree() for idx in range(10)]
        tree1 = trees[0]
        for tree2 in trees[1:]:
            x = test_target(
                    tree1=tree1,
                    tree2=tree2,
                    tree1_assemblage_leaf_sets=tree1.assemblage_leaf_sets,
                    tree2_assemblage_leaf_sets=tree2.assemblage_leaf_sets,
                    )
        self.validate_managed_and_cached_trees(
                test_target=test_target,
                trees=trees)

if __name__ == "__main__":
    unittest.main()

