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

from dendropy.datamodel import base

class AnnotableComparator(object):

    def compare_annotables(self, x1, x2):
        self.assertIsNot(x1, x2)
        if not x1.has_annotations:
            self.assertTrue( (not hasattr(x1, "_annotations")) or len(x1._annotations) == 0 )
            self.assertFalse(x2.has_annotations)
            self.assertTrue( (not hasattr(x2, "_annotations")) or len(x2._annotations) == 0 )
            return
        self.assertTrue( hasattr(x1, "_annotations") and len(x1._annotations) > 0 )
        self.assertTrue(x2.has_annotations)
        self.assertTrue( hasattr(x2, "_annotations") and len(x2._annotations) > 0 )
        self.assertIs(x1._annotations.target, x1)
        self.assertIs(x2._annotations.target, x2)
        self.assertIsNot(x1._annotations, x2._annotations)
        self.assertEqual(len(x1._annotations), len(x2._annotations))
        for a1, a2 in zip(x1._annotations, x2._annotations):
            self.assertIsNot(a1, a2)
            for k in a1.__dict__:
                self.assertIn(k, a2.__dict__)
                v1 = a1.__dict__[k]
                v2 = a2.__dict__[k]
                if isinstance(v1, base.DataObject):
                    self.assertTrue(isinstance(v2, base.DataObject))
                    self.assertIsNot(v1, v2)
                elif isinstance(v1, base.AnnotationSet):
                    self.assertTrue(isinstance(v2, base.AnnotationSet))
                    self.assertIs(v1.target, a1)
                    self.assertIs(v2.target, a2)
                    for s1, s2 in zip(v1, v2):
                        self.compare_annotables(s1, s2)
                else:
                    self.assertEqual(v1, v2)
                self.compare_annotables(a1, a2)
