import dendropy

# Create the DataSet to store data
ds = dendropy.DataSet()

# Set it up to manage all data under a single taxon namespace.
# HIGHLY RECOMMENDED!
taxon_namespace = dendropy.TaxonNamespace()
ds.attach_taxon_namespace(taxon_namespace)

# Read from multiple sources

# Add a collection of trees
ds.read(
    path='pythonidae.mle.nex',
    schema='nexus',)

# Add a collection of characters from a Nexus source
ds.read(
    path='pythonidae.chars.nexus',
    schema='nexus',)

# Add a collection of characters from a FASTA source
# Note that with this format, we have to explicitly provide the type of data
ds.read(
    path='pythonidae_cytb.fasta',
    schema='fasta',
    data_type="dna")

# Add a collection of characters from a PHYLIP source
# Note that with this format, we have to explicitly provide the type of data
ds.read(
    path='pythonidae.chars.phylip',
    schema='phylip',
    data_type="dna")

# Add a collection of continuous characters from a NeXML source
ds.read(
    path='pythonidae_continuous.chars.nexml',
    schema='nexml',)



