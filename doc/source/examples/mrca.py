import dendropy

tree = dendropy.Tree.get(path="pythonidae.mle.nex", schema="nexus")
taxon_labels=['Python sebae',
              'Python regius',
              'Python curtus',
              'Python molurus']
mrca = tree.mrca(taxon_labels=taxon_labels)
print(mrca.description())
