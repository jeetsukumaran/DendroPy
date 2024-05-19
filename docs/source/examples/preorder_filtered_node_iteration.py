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
multifurcating = lambda x: True if len(x.child_nodes()) > 2 else False
for nd in mle.postorder_node_iter(multifurcating):
    print(nd.description(0))
