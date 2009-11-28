#! /usr/bin/env python

import dendropy
from dendropy import treecalc
from dendropy import treesplit

tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")
taxon_labels=['Python sebae',
              'Python regius',
              'Python brongersmai',
              'Python molurus']
mrca = treecalc.find_mrca(tree, taxon_labels=taxon_labels)
print(mrca.description())
