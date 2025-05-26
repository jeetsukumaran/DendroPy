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
This subpackage handles the core definition of tree data structure class,
as well as all the structural classes that make up a tree.
"""

from dendropy.datamodel.treemodel._bipartition import Bipartition
from dendropy.datamodel.treemodel._edge import Edge
from dendropy.datamodel.treemodel._node import Node
from dendropy.datamodel.treemodel._tree import Tree
