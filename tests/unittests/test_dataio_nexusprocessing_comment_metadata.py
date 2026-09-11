#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Tests for comment metadata extraction (``dendropy.dataio.nexusprocessing``),
including support for a callable ``extract_comment_metadata`` argument and
the BEAST2-version-specific parser added to address
https://github.com/jeetsukumaran/DendroPy/issues/145: DendroPy's original
comment metadata parser mis-parses nested list-valued ("vector")
annotations, e.g. ``history_all={{57,0.08,C,T},{134,0.079,A,G}}``.
"""

import sys
import os
import unittest
import warnings
import dendropy
from dendropy.dataio import nexusprocessing

sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest

# Verbatim (whitespace preserved) from a node comment in a real BEAST/
# TreeAnnotator-format MCC tree (``HA_discrete_MCC.tre``, bundled as an
# example with the FigTree v1.4.4 source distribution).
REAL_MCC_TREE_EXAMPLE_COMMENT = (
        '&rate_range={4.938776751387227E-4,0.036916549293719556},'
        'height_median=9.000000000000002,'
        'length=0.6825227857160625,'
        'state="D",'
        'state.prob=1.0,'
        'rate=0.007334968720519001'
        )

# The nested-list annotation from issue #145.
ISSUE_145_COMMENT = "&history_all={{57,0.08,C,T},{134,0.079,A,G},{4,0.07,C,T}}"


def annotations_as_dict(annotations):
    return {a.name: a.value for a in annotations}


class DendroPyV5_0_0CommentMetadataParsingTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_dendropy_v5_0_0`` is DendroPy's original
    comment metadata parser, retained as the default: its nested-list
    limitation is part of the behavior tested for here.
    """

    def test_simple_scalar_values(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&rate=0.5,label=foo")
        self.assertEqual(d, [("rate", "0.5"), ("label", "foo")])

    def test_quoted_and_boolean_values(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                '&name="hello world",flag=true,off=FALSE')
        self.assertEqual(
                d, [("name", "hello world"), ("flag", True), ("off", False)])

    def test_single_level_vector(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("&x={1,2,3}")
        self.assertEqual(d, [("x", ["1", "2", "3"])])

    def test_nhx_format(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&&NHX:S=human:E=1.1.1.1")
        self.assertEqual(d, [("S", "human"), ("E", "1.1.1.1")])

    def test_nhx_format_with_explicit_prefix(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&&NHX:S=human")
        d2 = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("&&S=human")
        self.assertEqual(d, d2)

    def test_field_value_types_scalar(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&age=5", field_value_types={"age": float})
        self.assertEqual(d, [("age", 5.0)])

    def test_field_value_types_vector(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&x={1,2,3}", field_value_types={"x": int})
        self.assertEqual(d, [("x", [1, 2, 3])])

    def test_repeated_field_names_are_all_kept(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("&x=1,x=2,y=9")
        self.assertEqual(d, [("x", "1"), ("x", "2"), ("y", "9")])

    def test_unrecognized_comment_returns_empty(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("just a comment")
        self.assertEqual(d, [])

    def test_issue_145_nested_list_is_mis_parsed(self):
        # Regression-locks the originally-reported bug: the outer vector
        # value is truncated at the first (inner) closing brace.
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(ISSUE_145_COMMENT)
        self.assertEqual(d, [("history_all", ["{57", "0.08", "C", "T"])])


class ParseCommentMetadataToAnnotationsBackwardCompatTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_to_annotations`` retains its original
    signature and Annotation-set-returning behavior.
    """

    def test_basic(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations("&rate=0.5")
        self.assertEqual(annotations_as_dict(annotations), {"rate": "0.5"})

    def test_field_name_map(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                "&rate=0.5", field_name_map={"rate": "substitution_rate"})
        self.assertEqual(
                annotations_as_dict(annotations), {"substitution_rate": "0.5"})

    def test_field_value_types_on_vector(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                "&x={1,2,3}", field_value_types={"x": int})
        self.assertEqual(annotations_as_dict(annotations), {"x": [1, 2, 3]})

    def test_repeated_field_names_yield_one_annotation_each(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                "&x=1,x=2,y=9")
        self.assertEqual(
                sorted((a.name, a.value) for a in annotations),
                [("x", "1"), ("x", "2"), ("y", "9")])

    def test_accumulates_into_existing_set(self):
        existing = nexusprocessing.parse_comment_metadata_to_annotations("&a=1")
        combined = nexusprocessing.parse_comment_metadata_to_annotations(
                "&b=2", annotations=existing)
        self.assertIs(combined, existing)
        self.assertEqual(annotations_as_dict(combined), {"a": "1", "b": "2"})


class Beast2V2_7_8CommentMetadataParsingTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_beast2_v2_7_8`` reproduces the comment
    metadata parsing behavior of BEAST2 v2.7.8's ``TreeParser``,
    including correct handling of arbitrarily-nested list-valued
    annotations.
    """

    def test_numbers_and_strings(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                '&rate=0.0123,label="hello",tag=bareword')
        self.assertEqual(
                list(d), [("rate", 0.0123), ("label", "hello"), ("tag", "bareword")])

    def test_single_quoted_value(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&label='hello'")
        self.assertEqual(list(d), [("label", "hello")])

    def test_negative_and_scientific_notation_numbers(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                "&a=-1.5,b=4.938776751387227E-4")
        self.assertEqual(dict(d)["a"], -1.5)
        self.assertAlmostEqual(dict(d)["b"], 4.938776751387227E-4)

    def test_all_numeric_vector(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&hpd={1.1,2.2,3.3}")
        self.assertEqual(list(d), [("hpd", [1.1, 2.2, 3.3])])

    def test_mixed_vector_falls_back_to_raw_text(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8('&x={1,"a"}')
        self.assertEqual(list(d), [("x", ["1", '"a"'])])

    def test_string_only_vector_falls_back_to_raw_text(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&x={C,T,A,G}")
        self.assertEqual(list(d), [("x", ["C", "T", "A", "G"])])

    def test_nested_vector_resolves_issue_145(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(ISSUE_145_COMMENT)
        self.assertEqual(
                list(d),
                [("history_all",
                  ["{57,0.08,C,T}", "{134,0.079,A,G}", "{4,0.07,C,T}"])])

    def test_doubly_nested_vector(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&x={{{1,2},{3,4}},{5,6}}")
        self.assertEqual(list(d), [("x", ["{{1,2},{3,4}}", "{5,6}"])])

    def test_whitespace_is_tolerated(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                "& rate = 0.5 , hpd = { 1.1 , 2.2 } ")
        self.assertEqual(list(d), [("rate", 0.5), ("hpd", [1.1, 2.2])])

    def test_key_is_not_unquoted(self):
        # unlike a value, TreeParser.processMetadata() takes a key's
        # getText() as-is: a quoted key keeps its quotes (and so its
        # whitespace, which is otherwise insignificant)
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8('&" a key "=1')
        self.assertEqual(list(d), [('" a key "', 1.0)])
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&'a key'=1")
        self.assertEqual(list(d), [("'a key'", 1.0)])

    def test_bare_double_ampersand_is_not_stripped_and_does_not_warn(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                    "&&subject='Pythonidae'")
        self.assertEqual(caught, [])
        self.assertEqual(list(d), [("&subject", "Pythonidae")])

    def test_nhx_marker_warns(self):
        with self.assertWarns(UserWarning):
            d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                    "&&NHX:subject=Pythonidae")
        self.assertEqual(list(d), [("&NHX:subject", "Pythonidae")])

    def test_nhx_marker_warning_suggests_re_sub_workaround(self):
        import re
        extract_comment_metadata = lambda c: (
                nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                        re.sub(r"^&&NHX:?", "&", c)))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            extract_comment_metadata("&&NHX:subject=Pythonidae")
        self.assertEqual(caught, [])

    def test_nhx_marker_warning_suggests_whitespace_workaround(self):
        extract_comment_metadata = lambda c: (
                nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                        "& " + c[1:] if c.startswith("&&NHX") else c))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            extract_comment_metadata("&&NHX:subject=Pythonidae")
        self.assertEqual(caught, [])

    def test_empty_comment_returns_empty(self):
        self.assertEqual(
                list(nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&")), [])

    def test_repeated_field_names_keep_the_last_value(self):
        # BEAST2 calls node.setMetaData() per attribute, so a repeated
        # field name overwrites rather than accumulating
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&x=1,x=2,y=9")
        self.assertEqual(list(d), [("x", 2.0), ("y", 9.0)])

    def test_unrecognized_comment_returns_empty(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8("not a metadata comment")
        self.assertEqual(list(d), [])

    def test_malformed_missing_value_raises(self):
        with self.assertRaises(ValueError):
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&rate=")

    def test_malformed_unbalanced_vector_raises(self):
        with self.assertRaises(ValueError):
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&x={1,2")

    def test_malformed_trailing_content_raises(self):
        with self.assertRaises(ValueError):
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&rate=0.5,,")

    def test_numeric_only_key_raises(self):
        # per the BEAST2 grammar, an attribute key must lex as ASTRING,
        # not as a number
        with self.assertRaises(ValueError):
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8("&123=5")

    def test_real_mcc_tree_example_comment(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                REAL_MCC_TREE_EXAMPLE_COMMENT)
        d = dict(d)
        self.assertEqual(d["state"], "D")
        self.assertEqual(d["state.prob"], 1.0)
        self.assertAlmostEqual(d["rate"], 0.007334968720519001)
        self.assertEqual(
                d["rate_range"],
                [4.938776751387227E-4, 0.036916549293719556])


class Beast2V2_7_8NestingCommentMetadataParsingTestCase(dendropytest.ExtendedTestCase):

    def test_scalars_match_base_parser(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                '&rate=0.0123,label="hello",tag=bareword')
        self.assertEqual(
                list(d), [("rate", 0.0123), ("label", "hello"), ("tag", "bareword")])

    def test_flat_all_numeric_vector_matches_base_parser(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&hpd={1.1,2.2,3.3}")
        self.assertEqual(list(d), [("hpd", [1.1, 2.2, 3.3])])

    def test_flat_string_only_vector_matches_base_parser(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&x={C,T,A,G}")
        self.assertEqual(list(d), [("x", ["C", "T", "A", "G"])])

    def test_flat_mixed_vector_is_individually_typed(self):
        # the base parser falls back to raw text for the whole vector
        # here (["1", '"a"']); the nesting parser instead materializes
        # each element on its own
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                '&x={1,"a"}')
        self.assertEqual(list(d), [("x", [1.0, "a"])])

    def test_nested_numeric_vector_resolves_to_nested_floats(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&x={{1,2},{3,4}}")
        self.assertEqual(list(d), [("x", [[1.0, 2.0], [3.0, 4.0]])])

    def test_nested_vector_resolves_issue_145_with_nested_types(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                ISSUE_145_COMMENT)
        self.assertEqual(
                list(d),
                [("history_all",
                  [[57.0, 0.08, "C", "T"],
                   [134.0, 0.079, "A", "G"],
                   [4.0, 0.07, "C", "T"]])])

    def test_doubly_nested_all_numeric_vector_fully_resolves(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&x={{{1,2},{3,4}},{5,6}}")
        self.assertEqual(
                list(d), [("x", [[[1.0, 2.0], [3.0, 4.0]], [5.0, 6.0]])])

    def test_nested_vector_mixing_all_numeric_and_string_subvectors(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&x={{A,T},{57,0.08}}")
        self.assertEqual(list(d), [("x", [["A", "T"], [57.0, 0.08]])])

    def test_repeated_field_names_keep_the_last_value(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "&x=1,x=2,y=9")
        self.assertEqual(list(d), [("x", 2.0), ("y", 9.0)])

    def test_bare_double_ampersand_is_not_stripped_and_does_not_warn(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                    "&&subject='Pythonidae'")
        self.assertEqual(caught, [])
        self.assertEqual(list(d), [("&subject", "Pythonidae")])

    def test_nhx_marker_warns(self):
        with self.assertWarns(UserWarning):
            d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                    "&&NHX:subject=Pythonidae")
        self.assertEqual(list(d), [("&NHX:subject", "Pythonidae")])

    def test_empty_comment_returns_empty(self):
        self.assertEqual(
                list(nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting("&")),
                [])

    def test_unrecognized_comment_returns_empty(self):
        d = nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                "not a metadata comment")
        self.assertEqual(list(d), [])

    def test_malformed_unbalanced_vector_raises(self):
        with self.assertRaises(ValueError):
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting("&x={1,2")

    def test_real_mcc_tree_example_comment_has_no_nested_vectors(self):
        # this fixture has no nested vectors, so both parsers agree
        base = dict(nexusprocessing.parse_comment_metadata_beast2_v2_7_8(
                REAL_MCC_TREE_EXAMPLE_COMMENT))
        nested = dict(nexusprocessing.parse_comment_metadata_beast2_v2_7_8_nesting(
                REAL_MCC_TREE_EXAMPLE_COMMENT))
        self.assertEqual(base, nested)


class CommentMetadataToAnnotationsTestCase(dendropytest.ExtendedTestCase):

    def test_basic(self):
        annotations = nexusprocessing._comment_metadata_to_annotations(
                [("rate", 0.5), ("label", "x")])
        self.assertEqual(
                annotations_as_dict(annotations), {"rate": 0.5, "label": "x"})

    def test_field_name_map(self):
        annotations = nexusprocessing._comment_metadata_to_annotations(
                [("rate", 0.5)], field_name_map={"rate": "substitution_rate"})
        self.assertEqual(
                annotations_as_dict(annotations), {"substitution_rate": 0.5})

    def test_accepts_pairs_with_repeated_field_names(self):
        annotations = nexusprocessing._comment_metadata_to_annotations(
                [("x", 1), ("x", 2)])
        self.assertEqual(
                sorted((a.name, a.value) for a in annotations),
                [("x", 1), ("x", 2)])

    def test_empty_metadata_yields_no_annotations(self):
        annotations = nexusprocessing._comment_metadata_to_annotations([])
        self.assertEqual(len(annotations), 0)


class ExtractCommentMetadataCallableIntegrationTestCase(dendropytest.ExtendedTestCase):
    """
    End-to-end tests confirming that a callable ``extract_comment_metadata``
    is honored by the tree readers (for both NEWICK and NEXUS schemas),
    and that using the BEAST2-style parser resolves issue #145.
    """

    NEWICK_STR = (
            "((A[&rate=0.5,hpd={1.1,2.2}]:1,"
            "B[&history_all={{57,0.08,C,T},{134,0.079,A,G}}]:1):1,C:1);"
            )

    def _annotations_by_taxon_label(self, tree):
        result = {}
        for nd in tree:
            if nd.taxon is not None:
                result[nd.taxon.label] = nd.annotations.values_as_dict()
        return result

    def test_default_bool_true_matches_dendropy_v5_0_0(self):
        tree = dendropy.Tree.get(data=self.NEWICK_STR, schema="newick")
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(result["A"]["rate"], "0.5")
        self.assertEqual(result["B"]["history_all"], ["{57", "0.08", "C", "T"])

    def test_extract_comment_metadata_false(self):
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=False)
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "B":
                self.assertEqual(len(nd.comments), 1)
                self.assertIn("history_all", nd.comments[0])
                self.assertEqual(len(nd.annotations), 0)

    def test_extract_comment_metadata_callable_beast2_resolves_issue_145(self):
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=nexusprocessing.parse_comment_metadata_beast2_v2_7_8)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(
                result["B"]["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}"])
        self.assertEqual(result["A"]["rate"], 0.5)
        self.assertEqual(result["A"]["hpd"], [1.1, 2.2])

    def test_extract_comment_metadata_callable_via_nexus_schema(self):
        nexus_str = "#NEXUS\nBegin trees;\n  tree t1 = " + self.NEWICK_STR + "\nEnd;\n"
        tree = dendropy.Tree.get(
                data=nexus_str, schema="nexus",
                extract_comment_metadata=nexusprocessing.parse_comment_metadata_beast2_v2_7_8)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(
                result["B"]["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}"])

    def test_annotations_follow_comment_order(self):
        # annotation order is stable from run to run: the parsed pairs
        # are applied to the target set directly, rather than by way of
        # an intermediate (unordered) ``set``
        tree = dendropy.Tree.get(
                data="(A[&aa=1,bb=2,cc=3,dd=4,ee=5,ff=6]:1,B:1);",
                schema="newick")
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "A":
                self.assertEqual(
                        [a.name for a in nd.annotations],
                        ["aa", "bb", "cc", "dd", "ee", "ff"])

    def test_custom_user_callable(self):
        def my_parser(comment):
            return [("raw", comment)]
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=my_parser)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(result["A"]["raw"], "&rate=0.5,hpd={1.1,2.2}")

    def test_custom_callable_returning_generator(self):
        # a generator is always truthy, so the emptiness test must look
        # at the pairs themselves rather than the iterable
        def empty_gen(comment):
            return (x for x in [])
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=empty_gen)
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "A":
                self.assertEqual(len(nd.annotations), 0)
                self.assertEqual(len(nd.comments), 1)

    def test_custom_callable_returning_empty_falls_back_to_comments(self):
        # when metadata extraction yields nothing, the raw comment is
        # kept (as a plain comment) rather than silently dropped
        def no_op_parser(comment):
            return []
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=no_op_parser)
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "A":
                self.assertEqual(len(nd.annotations), 0)
                self.assertEqual(len(nd.comments), 1)


if __name__ == "__main__":
    unittest.main()
