#! /usr/bin/env python

import dendropy

tree = dendropy.Tree.get_from_path("pythonidae.mle.nex", "nexus")
taxon_labels=['Python sebae',
              'Python regius',
              'Python curtus',
              'Python molurus']
mrca = tree.mrca(taxon_labels=taxon_labels)
print(mrca.description())
