#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import dendropy

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)

mle = dendropy.Tree.get(
    path='pythonidae.mle.nex',
    schema='nexus')
short = lambda edge: True if edge.length < 0.01 else False
for edge in mle.postorder_edge_iter(short):
    print(edge.length)
