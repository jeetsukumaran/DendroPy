import warnings
import dendropy

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)


ds1 = dendropy.DataSet()
taxon_namespace = dendropy.TaxonNamespace()
ds1.attach_taxon_namespace(taxon_namespace)
ds1.read(
    path='pythonidae.mle.nex',
    schema='nexus',)
ds1.read(
    path='pythonidae.chars.nexus',
    schema='nexus',)

ds2 = dendropy.DataSet(ds1)

# for tx1, tx2 in zip(ds1, ds2):
#     print(tx1 is tx2)
