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
Wrapper for interacting with RSPR
"""

from io import StringIO
import subprocess
import socket
import os

from dendropy.utility.messaging import get_logger
from dendropy.utility import processio
_LOG = get_logger("interop.rspr")

HOSTNAME = socket.gethostname()
PID = os.getpid()

class Rspr(object):
    """
    This class wraps all attributes and input needed to make a call to RSPR.

        https://github.com/cwhidden/rspr

    RSPR:

    Calculate approximate and exact Subtree Prune and Regraft (rSPR)
    distances and the associated maximum agreement forests (MAFs) between pairs
    of rooted binary trees from STDIN in newick format. Supports arbitrary labels.
    The second tree may be multifurcating.

    Copyright 2009-2014 Chris Whidden
    whidden@cs.dal.ca
    http://kiwi.cs.dal.ca/Software/RSPR

    """
    ###
    # NOTE ON THE ``--pairwise`` FLAG
    # -------------------------------
    # This determines the type of comparisons done.
    #
    # Specified without arguments, it does all distinct pairwise comparisons of
    # the input tree. Output by default is a matrix with only the upper half
    # filled. So, assuming the source has 4 trees, then:
    #
    #   $ cat 4trees.tre | rspr -pairwise
    #   0,23,24,24
    #   ,0,5,7
    #   ,,0,6
    #
    # Further (numerical) arguments specify the row/columns to restrict the
    # comparisons.
    #
    # E.g., first row only: first tree to all other trees:
    #
    #   $ cat 4trees.tre | rspr -pairwise 0 1
    #   0,23,24,24
    #
    # E.g. first two rows only: first two trees to all other trees:
    #
    #   $ cat 4trees.tre | rspr -pairwise 0 2
    #   0,23,24,24
    #   ,0,5,7
    #
    # E.g. third row only:
    #
    #   $ cat 4trees.tre | rspr -pairwise 2 3
    #   ,,0,6
    #
    # E.g. last column only: last tree to all other trees:
    #
    #   $ cat 4trees.tre | rspr -pairwise 0 4 3 4
    #   24
    #   7
    #   6
    #   0


    def __init__(self,
            algorithm="bb",
            optimizations=None,
            cc=None,
            ):
        """

        Parameters
        ----------
        algorithm : str
            One of "fpt", "bb", "approx".
        optimizations: list[str] or None
            Will be passed directly rspr (with each element prefixed by ``-``).
        cc: bool
            Calculate a potentially better approximation with a quadratic time
            algorithm.
        """
        self.algorithm = algorithm
        self.optimizations = optimizations

    def compare_one_to_many(self,
            ref_tree,
            comparison_trees,
            command_args=None,
            newick_output_kwargs=None,
            ):
        """

        Compare ``ref_tree'' to each tree in ``comparison_trees``.

        Parameters
        ----------
        ref_tree : |Tree|
            A |Tree| object to be compared to every tree in ``comparison_trees``.
        comparison_trees : |Tree|
            An (ordered) iterable of trees to which ``ref_tree`` should be
            compared.
        command_args : list or None
            An iterable of (string) arguments to be passed to the program.
        newick_output_kwargs : dict or None
            A collection of keyword arguments to pass to the tree string
            composition routines (that will generate the tree strings to be
            used as input to rspr).

        Returns
        -------
        scores : list[numeric]
            A list of the SPR distances from ``ref_tree'' to
            ``comparison_trees``, in order of the trees given.
        """
        if newick_output_kwargs is None:
            newick_output_kwargs = {}
        tf = StringIO()
        ref_tree.write(file=tf, schema="newick", **newick_output_kwargs)
        for t in comparison_trees:
            t.write(file=tf, schema="newick", **newick_output_kwargs)
        command = []
        command.append("rspr") # TODO: command path as instance attribute
        command.extend(["-pairwise", "0", "1"])
        if command_args is not None:
            command.extend(command_args)
        p = subprocess.Popen(command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,)
        stdout, stderr = processio.communicate(p, commands=tf.getvalue())
        result_fields = stdout.strip("\n").split(",")
        assert len(result_fields) == 1 + len(comparison_trees), "Expecting length {} + 1 for results, but received {}: {}".format(len(comparison_trees), len(result_fields), result_fields)
        return [int(v) for v in result_fields[1:]]








