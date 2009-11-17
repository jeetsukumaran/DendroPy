*************************
Taxa and Taxon Management
*************************

Operational taxonomic units in DendroPy are represented by :class:`~dendropy.dataobject.taxon.Taxon` objects, and distinct collections of operational taxonomic units are represented by :class:`~dendropy.dataobject.taxon.TaxonSet` objects.
Every time a definition of taxa is encountered in a data source, for example, a "TAXA" block in a NEXUS file, a new :class:`~dendropy.dataobject.taxon.TaxonSet` object is created and populated with :class:`~dendropy.dataobject.taxon.Taxon` objects corresponding to the taxa defined in the data source.
Some data formats do not have explicit definition of taxa, e.g. a NEWICK tree source.
These nonetheless can be considered to have an implicit definition of a collection of operational taxonomic units given by the aggregate of all operational taxonomic units referenced in the data (i.e., the set of all distinct labels on trees in a NEWICK file).
