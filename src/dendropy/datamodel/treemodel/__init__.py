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
This subpackage handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from dendropy.datamodel.treemodel._bipartition import Bipartition
from dendropy.datamodel.treemodel._edge import Edge
from dendropy.datamodel.treemodel._node import Node
from dendropy.datamodel.treemodel._tree import Tree
