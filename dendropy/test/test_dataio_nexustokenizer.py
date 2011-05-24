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
NEXUS data read/write parse/format tests.
"""

import sys
import os
import unittest
from cStringIO import StringIO

from dendropy.utility import messaging
from dendropy.dataio import nexustokenizer
import dendropy

_LOG = messaging.get_logger(__name__)

class NexusTokenTest(unittest.TestCase):

    def testParseSpacy(self):
        f = """#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=30;
    TAXLABELS
        'Anolis ahli'
        'Anolis garmani'
        'Anolis grahami'
        'Anolis valencienni'
        'Anolis lineatopus'
        'Anolis aliniger'
        'Anolis coelestinus'
        'Anolis bahorucoensis'
        'Anolis equestris'
        'Anolis luteogularis'
        'Anolis occultus'
        'Anolis barahonae'
        'Anolis cuvieri'
        'Anolis insolitus'
        'Anolis olssoni'
        'Anolis brevirostris'
        'Anolis distichus'
        'Anolis cristatellus'
        'Anolis krugi'
        'Anolis stratulus'
        'Anolis alutaceus'
        'Anolis vanidicus'
        'Anolis angusticeps'
        'Anolis paternus'
        'Anolis loysiana'
        'Anolis marcanoi'
        'Anolis strahmi'
        'Diplolaemus darwinii'
        'Anolis ophiolepis'
        'Anolis sagrei'
  ;
END;

BEGIN TREES;
    tree 'con 50 majrule' = [&U] ('Anolis ahli':0.2642130000,((('Anolis garmani':0.1068380000,'Anolis grahami':0.0863670000)1.00:0.069511,'Anolis valencienni':0.1642630000)0.87:0.020752,'Anolis lineatopus':0.1957260000)1.00:0.077682,((((((('Anolis aliniger':0.1600010000,'Anolis coelestinus':0.1932310000)1.00:0.071920,'Anolis bahorucoensis':0.2266880000)0.68:0.023043,('Anolis equestris':0.0227020000,'Anolis luteogularis':0.0306410000)1.00:0.198165,'Anolis occultus':0.4231200000)0.89:0.056277,('Anolis barahonae':0.2114890000,'Anolis cuvieri':0.1686700000)1.00:0.084190,('Anolis insolitus':0.2438820000,'Anolis olssoni':0.2568770000)1.00:0.050618)0.86:0.031679,(('Anolis brevirostris':0.1801300000,'Anolis distichus':0.1151360000)1.00:0.123136,(('Anolis cristatellus':0.2144360000,'Anolis krugi':0.1573300000)0.93:0.036788,'Anolis stratulus':0.1973470000)1.00:0.081037)1.00:0.056582)0.77:0.021826,(('Anolis alutaceus':0.1619060000,'Anolis vanidicus':0.2059960000)1.00:0.118216,(('Anolis angusticeps':0.0857100000,'Anolis paternus':0.0595110000)1.00:0.153413,'Anolis loysiana':0.1836280000)1.00:0.042858)1.00:0.057139,('Anolis marcanoi':0.2359120000,'Anolis strahmi':0.1977660000)1.00:0.141032,'Diplolaemus darwinii':0.6364930000)1.00:0.067869,('Anolis ophiolepis':0.0945010000,'Anolis sagrei':0.0967580000)1.00:0.179398)0.96:0.044895);
END;
"""
        r2 = StringIO(f)
        tokenizer = nexustokenizer.NexusTokenizer(r2)
        expected = ['#NEXUS', 'BEGIN', 'TAXA', ';', 'DIMENSIONS', 'NTAX', '=',
                    '30', ';', 'TAXLABELS', 'Anolis ahli', 'Anolis garmani',
                    'Anolis grahami', 'Anolis valencienni', 'Anolis lineatopus',
                    'Anolis aliniger', 'Anolis coelestinus', 'Anolis bahorucoensis',
                    'Anolis equestris', 'Anolis luteogularis', 'Anolis occultus', 'Anolis barahonae',
                    'Anolis cuvieri', 'Anolis insolitus', 'Anolis olssoni', 'Anolis brevirostris',
                    'Anolis distichus', 'Anolis cristatellus', 'Anolis krugi', 'Anolis stratulus', 'Anolis alutaceus', 'Anolis vanidicus', 'Anolis angusticeps', 'Anolis paternus', 'Anolis loysiana', 'Anolis marcanoi', 'Anolis strahmi', 'Diplolaemus darwinii', 'Anolis ophiolepis', 'Anolis sagrei',
                    ';', 'END', ';', 'BEGIN', 'TREES', ';', 'tree', 'con 50 majrule',  '=', '(', 'Anolis ahli', ':', '0.2642130000', ',', '(', '(', '(', 'Anolis garmani', ':', '0.1068380000', ',', 'Anolis grahami', ':', '0.0863670000', ')', '1.00', ':', '0.069511', ',', 'Anolis valencienni', ':', '0.1642630000', ')', '0.87', ':', '0.020752', ',', 'Anolis lineatopus', ':', '0.1957260000', ')', '1.00', ':', '0.077682', ',', '(', '(', '(', '(', '(', '(', '(', 'Anolis aliniger', ':', '0.1600010000', ',', 'Anolis coelestinus', ':', '0.1932310000', ')', '1.00', ':', '0.071920', ',', 'Anolis bahorucoensis', ':', '0.2266880000', ')', '0.68', ':', '0.023043', ',', '(', 'Anolis equestris', ':', '0.0227020000', ',', 'Anolis luteogularis', ':', '0.0306410000', ')', '1.00', ':', '0.198165', ',', 'Anolis occultus', ':', '0.4231200000', ')', '0.89', ':', '0.056277', ',', '(', 'Anolis barahonae', ':', '0.2114890000', ',', 'Anolis cuvieri', ':', '0.1686700000', ')', '1.00', ':', '0.084190', ',', '(', 'Anolis insolitus', ':', '0.2438820000', ',', 'Anolis olssoni', ':', '0.2568770000', ')', '1.00', ':', '0.050618', ')', '0.86', ':', '0.031679', ',', '(', '(', 'Anolis brevirostris', ':', '0.1801300000', ',', 'Anolis distichus', ':', '0.1151360000', ')', '1.00', ':', '0.123136', ',', '(', '(', 'Anolis cristatellus', ':', '0.2144360000', ',', 'Anolis krugi', ':', '0.1573300000', ')', '0.93', ':', '0.036788', ',', 'Anolis stratulus', ':', '0.1973470000', ')', '1.00', ':', '0.081037', ')', '1.00', ':', '0.056582', ')', '0.77', ':', '0.021826', ',', '(', '(', 'Anolis alutaceus', ':', '0.1619060000', ',', 'Anolis vanidicus', ':', '0.2059960000', ')', '1.00', ':', '0.118216', ',', '(', '(', 'Anolis angusticeps', ':', '0.0857100000', ',', 'Anolis paternus', ':', '0.0595110000', ')', '1.00', ':', '0.153413', ',', 'Anolis loysiana', ':', '0.1836280000', ')', '1.00', ':', '0.042858', ')', '1.00', ':', '0.057139', ',', '(', 'Anolis marcanoi', ':', '0.2359120000', ',', 'Anolis strahmi', ':', '0.1977660000', ')', '1.00', ':', '0.141032', ',', 'Diplolaemus darwinii', ':', '0.6364930000', ')', '1.00', ':', '0.067869', ',', '(', 'Anolis ophiolepis', ':', '0.0945010000', ',', 'Anolis sagrei', ':', '0.0967580000', ')', '1.00', ':', '0.179398', ')', '0.96', ':', '0.044895', ')', ';', 'END',';']
        for e in expected:
            token = tokenizer.read_next_token()
            self.assertEqual(token, e)

    def testComments(self):
        f = "  i [f]t [i] y"
        r2 = StringIO(f)
        tokenizer = nexustokenizer.NexusTokenizer(r2)
        expected = ['i', 't', 'y']
        for e in expected:
            token = tokenizer.read_next_token()
            self.assertEqual(e, token)

class CommentReadingTests(unittest.TestCase):

    def testSimplePostNodeComments(self):
        s = "((A[A]:1,B[B]:1)AB[AB]:1,(C[C]:1,D[D]:1)CD[CD]:1)Root[Root]:1"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.assign_node_labels_from_taxon_or_oid()
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.edge.length, 1)

    def testSimplePostEdgeLengthComments(self):
        s = "((A:1[A],B:1[B])AB:1[AB],(C:1[C],D:1[D])CD:1[CD])Root:1[Root]"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.assign_node_labels_from_taxon_or_oid()
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)

    def testPostNodeAndEdgeLengthComments(self):
        s = "((A[A]:1[A],B[B]:1[B])AB[AB]:1[AB],(C[C]:1[C],D[D]:1[D])CD[CD]:1[CD])Root[Root]:1[Root]"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.assign_node_labels_from_taxon_or_oid()
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 2)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.comments[1], nd.label)

    def testMultiPositionComments(self):
        s = """(([xxx]A[A][A]:[A][A]1[A][A],
                 [xxx]B[B][B]:[B][B]1[B][B])
                 [xxx]AB[AB][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[C][C]:[C][C]1[C][C],
                 [xxx]D[D][D]:[D][D]1[D][D])
                 [xxx]CD[CD][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[Root][Root]:[Root][Root]1[Root][Root]"""
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.assign_node_labels_from_taxon_or_oid()
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 7)
            self.assertEqual(nd.comments[0], 'xxx')
            for i in xrange(1,7):
                self.assertEqual(nd.comments[i], nd.label)

    def testIncompleteMetadata(self):
        s = """[&color=blue](A[&region=Asia,id=00012][cryptic],(B[&region=Africa],C[&region=Madagascar,id=19391][two of three]))"""
        tree = dendropy.Tree.get_from_string(s, 'newick', extract_comment_metadata=True)
        self.assertEqual(tree.comment_metadata, {'color': 'blue'})
        expected = [ {'region': 'Asia', 'id': '00012'},
                {'region': 'Africa'},
                {'region': 'Madagascar', 'id': '19391'},
                {},
                {},]
        for idx, nd in enumerate(tree.postorder_node_iter()):
            self.assertEqual(nd.comment_metadata, expected[idx])

class CommentMetaDataTests(unittest.TestCase):

    nexus_skeleton = """\
#NEXUS

BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS
        A
        B
        C
        D
  ;
END;

BEGIN TREES;
    TREE 0 = %s ;
END;
"""
    figtree_metadata_str = """[&Tree1=1,Tree2=2, Tree3={1,2,3}](([xxx]A[&A1=1,A2=2,A3={1,2,3},  ,][A]:[A][A]1[A][A],
                 [xxx]B[&B1=1,B2=2,B3={1,2,3}][B]:[B][B]1[B][B])
                 [xxx]AB[&AB1=1,AB2=2,AB3={1,2,3}][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[&C1=1,C2=2,C3={1,2,3}][C]:[C][C]1[C][C],
                 [xxx]D[&D1=1,D2=2,D3={1,2,3}][D]:[D][D]1[D][D])
                 [xxx]CD[&CD1=1, CD2=2, CD3={1,2,3}][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[&Root1=1, Root2=2, Root3={1,2,3}][Root]:[Root][Root]1[Root][Root]"""

    nhx_metadata_str = """[&Tree1=1,Tree2=2, Tree3={1,2,3}](([xxx]A[&&A1=1:A2=2:A3={1,2,3}][A]:[A][A]1[A][A],
                 [xxx]B[&&B1=1:B2=2:B3={1,2,3}][B]:[B][B]1[B][B])
                 [xxx]AB[&&AB1=1:AB2=2:AB3={1,2,3}][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[&&C1=1:C2=2:C3={1,2,3}][C]:[C][C]1[C][C],
                 [xxx]D[&&D1=1:D2=2:D3={1,2,3}][D]:[D][D]1[D][D])
                 [xxx]CD[&&CD1=1: CD2=2: CD3={1,2,3}][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[&&Root1=1: Root2=2: Root3={1,2,3}][Root]:[Root][Root]1[Root][Root]"""

    def check_results(self, tree):
        tree.assign_node_labels_from_taxon_or_oid()
        self.assertEqual(tree.comment_metadata, {'Tree1': '1', 'Tree2': '2', 'Tree3':['1','2','3']})
        for nd in tree.postorder_node_iter():
            metadata = nd.comment_metadata
            #print("%s: %s => %s" % (nd.label, nd.comments, metadata))
            self.assertEqual(len(metadata), 3)
            values = ["1", "2", ["1","2","3"]]
            for i in range(3):
                key = "%s%d" % (nd.label, i+1)
                self.assertTrue(key in metadata)
                self.assertEqual(metadata[key], values[i])

    def testFigtreeStyleBasic(self):
        s = self.figtree_metadata_str
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.comment_metadata = nexustokenizer.parse_comment_metadata(tree.comments)
        for nd in tree.postorder_node_iter():
            nd.comment_metadata = nexustokenizer.parse_comment_metadata(nd.comments)
        self.check_results(tree)

    def testNHXBasic(self):
        s = self.nhx_metadata_str
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(s, 'newick')
        tree.comment_metadata = nexustokenizer.parse_comment_metadata(tree.comments)
        for nd in tree.postorder_node_iter():
            nd.comment_metadata = nexustokenizer.parse_comment_metadata(nd.comments)
        self.check_results(tree)

    def testNHXStyleNexusReader(self):
        s = self.nexus_skeleton % self.nhx_metadata_str
        tree = dendropy.Tree.get_from_string(s,
                'nexus',
                extract_comment_metadata=True)
        self.check_results(tree)

    def testFigtreeStyleNexusReader(self):
        s = self.nexus_skeleton % self.figtree_metadata_str
        tree = dendropy.Tree.get_from_string(s,
                'nexus',
                extract_comment_metadata=True)
        self.check_results(tree)

if __name__ == "__main__":
    unittest.main()
