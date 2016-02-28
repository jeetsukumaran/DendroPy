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

from . import pathmap
from dendropy.datamodel import basemodel
from dendropy.datamodel import charmatrixmodel
from dendropy.utility.textprocessing import StringIO

class ValidateWriteable(object):

    def write_out_validate_equal_and_return(self,
            writeable,
            schema,
            kwargs):
        self.maxDiff = None
        sio = StringIO()
        writeable.write(file=sio, schema=schema, **kwargs)
        s0 = sio.getvalue()
        s1 = writeable.as_string(schema=schema, **kwargs)
        self.assertEqual(s0, s1)
        with pathmap.SandboxedFile() as tempf:
            writeable.write(path=tempf.name, schema=schema, **kwargs)
            tempf.flush()
            tempf.close()
            with open(tempf.name, "r") as src:
                s2 = src.read()
        self.assertEqual(s0, s2)
        return s0

class ValidateReadable(object):

    def get_from(self,
            object_type,
            source_path,
            schema,
            kwargs):
        results = []
        with open(source_path, "r") as s0:
            obj0 = object_type.get(stream=s0, schema=schema, **kwargs)
            results.append(obj0)
        obj1 = object_type.get(path=source_path, schema=schema, **kwargs)
        results.append(obj1)
        with open(source_path, "r") as s0:
            s = s0.read()
            obj2 = object_type.get(data=s, schema=schema, **kwargs)
        results.append(obj2)
        return results

    def read_from(self,
            object_factory,
            source_path,
            schema,
            kwargs):
        results = []
        obj0 = object_factory()
        with open(source_path, "r") as s0:
            obj0.read(stream=s0, schema=schema, **kwargs)
            results.append(obj0)
        obj1 = object_factory()
        obj1.read(path=source_path, schema=schema, **kwargs)
        results.append(obj1)
        obj2 = object_factory()
        with open(source_path, "r") as s0:
            s = s0.read()
            obj2.read(value=s, schema=schema, **kwargs)
        results.append(obj2)
        return results

class Comparator(object):

    def compare_distinct_taxon(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_annotations=True):
        if taxon_namespace_scoped:
            self.assertIs(x1, x2)
        else:
            if x1 is None or x2 is None:
                self.assertIs(x1, None)
                self.assertIs(x2, None)
            else:
                self.assertIsNot(x1, x2)
                self.assertEqual(x1.label, x2.label)
                if compare_annotations:
                    self.compare_distinct_annotables(x1, x2)

    def compare_distinct_taxon_namespace(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_annotations=True,):
        if taxon_namespace_scoped:
            self.assertIs(x1, x2)
        else:
            if x1 is None or x2 is None:
                self.assertIs(x1, None)
                self.assertIs(x2, None)
            else:
                self.assertIsNot(x1, x2)
                self.assertEqual(x1.label, x2.label)
                self.assertEqual(len(x1._taxa), len(x2._taxa))
                for t1, t2 in zip(x1._taxa, x2._taxa):
                    self.compare_distinct_taxon(t1, t2,
                            taxon_namespace_scoped=taxon_namespace_scoped,
                            compare_annotations=compare_annotations)
                if compare_annotations:
                    self.compare_distinct_annotables(x1, x2)

    def compare_distinct_sequences(self,
            x1, x2,
            compare_annotations=True):
        self.assertIsNot(x1, x2)
        self.assertEqual(len(x1), len(x2))
        for c1, c2 in zip(x1, x2):
            self.assertEqual(c1, c2)
        if compare_annotations:
            self.compare_distinct_annotables(x1, x2)

    def compare_distinct_char_matrix(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_matrix_annotations=True,
            compare_sequence_annotations=True,
            compare_taxon_annotations=True):
        self.assertIsNot(x1, x2)
        self.assertEqual(len(x1), len(x2))
        self.compare_distinct_taxon_namespace(x1.taxon_namespace, x2.taxon_namespace,
                taxon_namespace_scoped=taxon_namespace_scoped,
                compare_annotations=compare_taxon_annotations,
                )
        self.assertEqual(x1.label, x2.label)
        self.assertEqual(x1.data_type, x2.data_type)
        if isinstance(x1, charmatrixmodel.DiscreteCharacterMatrix):
            # print(x1.data_type, x1.datatype_alphabet)
            self.assertEqual(len(x1.state_alphabets), len(x2.state_alphabets))
            for sa1, sa2 in zip(x1.state_alphabets, x2.state_alphabets):
                self.assertIs(sa1, sa2)
            self.assertIs(x1._default_state_alphabet, x2._default_state_alphabet)
        for t1, t2 in zip(x1, x2):
            if taxon_namespace_scoped:
                self.assertIs(t1, t2)
            else:
                self.assertIsNot(t1, t2)
            s1 = x1[t1]
            s2 = x2[t2]
            self.compare_distinct_sequences(s1, s2,
                    compare_annotations=compare_sequence_annotations)
        if compare_matrix_annotations:
            self.compare_distinct_annotables(x1, x2)

    def compare_distinct_tree_list(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_tree_annotations=True,
            compare_taxon_annotations=True):
        self.assertIsNot(x1, x2)
        self.assertEqual(len(x1), len(x2))
        self.compare_distinct_taxon_namespace(x1.taxon_namespace, x2.taxon_namespace,
                taxon_namespace_scoped=taxon_namespace_scoped,
                compare_annotations=compare_taxon_annotations,
                )
        self.assertEqual(x1.label, x2.label)
        for t1, t2 in zip(x1, x2):
            self.assertIsNot(t1, t2)
            self.compare_distinct_trees(t1, t2,
                    taxon_namespace_scoped=taxon_namespace_scoped,
                    compare_tree_annotations=compare_tree_annotations,
                    compare_taxon_annotations=compare_taxon_annotations)

    def compare_distinct_trees(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_tree_annotations=True,
            compare_taxon_annotations=True):
        self.assertIsNot(x1, x2)
        self.assertEqual(x1.label, x2.label)
        self.compare_distinct_taxon_namespace(x1.taxon_namespace, x2.taxon_namespace,
                taxon_namespace_scoped=taxon_namespace_scoped,
                compare_annotations=compare_taxon_annotations,
                )
        self.compare_distinct_nodes(x1.seed_node, x2.seed_node,
                taxon_namespace_scoped=taxon_namespace_scoped,
                compare_tree_annotations=compare_tree_annotations,
                compare_taxon_annotations=compare_taxon_annotations)

    def compare_distinct_nodes(self,
            x1, x2,
            taxon_namespace_scoped=True,
            compare_tree_annotations=True,
            compare_taxon_annotations=True,
            check_children=True):
        self.assertIsNot(x1, x2)
        self.assertEqual(x1.label, x2.label)
        taxon1 = x1.taxon
        taxon2 = x2.taxon
        self.compare_distinct_taxon(x1.taxon, x2.taxon,
                taxon_namespace_scoped=taxon_namespace_scoped,
                compare_annotations=compare_taxon_annotations)
        self.assertIsNot(x1.edge, x2.edge)
        self.assertEqual(x1.edge.label, x2.edge.label)
        self.assertEqual(x1.edge.length, x2.edge.length)
        self.assertIs(x1.edge.head_node, x1)
        self.assertIs(x2.edge.head_node, x2)
        if x1._parent_node is None or x2._parent_node is None:
            self.assertIs(x1._parent_node, None)
            self.assertIs(x2._parent_node, None)
            self.assertIs(x1.edge.tail_node, None)
            self.assertIs(x2.edge.tail_node, None)
        else:
            self.compare_distinct_nodes(x1._parent_node, x2._parent_node,
                    taxon_namespace_scoped=taxon_namespace_scoped,
                    compare_tree_annotations=compare_tree_annotations,
                    compare_taxon_annotations=compare_taxon_annotations,
                    check_children=False)
            self.assertIn(x1, x1._parent_node._child_nodes)
            self.assertNotIn(x1, x2._parent_node._child_nodes)
            self.assertIn(x2, x2._parent_node._child_nodes)
            self.assertNotIn(x2, x1._parent_node._child_nodes)
        if check_children:
            self.assertEqual(len(x1._child_nodes), len(x2._child_nodes))
            for c1, c2 in zip(x1._child_nodes, x2._child_nodes):
                self.compare_distinct_nodes(c1, c2,
                        taxon_namespace_scoped=taxon_namespace_scoped,
                        compare_tree_annotations=compare_tree_annotations,
                        compare_taxon_annotations=compare_taxon_annotations,
                        check_children=True)
        if compare_tree_annotations:
            self.compare_distinct_annotables(x1, x2)
            self.compare_distinct_annotables(x1.edge, x2.edge)

    def compare_distinct_annotables(self, x1, x2):
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
            if a1.is_attribute:
                self.assertTrue(a2.is_attribute)
                self.assertEqual(a1._value[1], a2._value[1])
            else:
                self.assertEqual(a1._value, a2._value)
            for k in a1.__dict__:
                if k == "_value":
                    continue
                self.assertIn(k, a2.__dict__)
                v1 = a1.__dict__[k]
                v2 = a2.__dict__[k]
                if isinstance(v1, basemodel.DataObject):
                    self.assertTrue(isinstance(v2, basemodel.DataObject))
                    self.assertIsNot(v1, v2)
                elif isinstance(v1, basemodel.AnnotationSet):
                    self.assertTrue(isinstance(v2, basemodel.AnnotationSet))
                    self.assertIs(v1.target, a1)
                    self.assertIs(v2.target, a2)
                    for s1, s2 in zip(v1, v2):
                        self.compare_distinct_annotables(s1, s2)
                else:
                    self.assertEqual(v1, v2)
                self.compare_distinct_annotables(a1, a2)
