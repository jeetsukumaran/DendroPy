*******************************
Converting Between Data Formats
*******************************

Any data in a schema that can be read by DendroPy, can be saved to files in any schema that can be written by DendroPy.
Converting data between formats is simply a matter of calling readers and writers of the appropriate type.

Converting from FASTA schema to NEXUS::

    >>> import dendropy
    >>> cytb = dendropy.DnaCharacterMatrix.get_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> cytb.write_to_path("pythonidae_cytb.nexus", "nexus")

Converting a collection of trees from NEXUS schema to Newick::

    >>> import dendropy
    >>> mcmc = dendropy.TreeList.get_from_path("pythonidae.mcmc.nex", "nexus")
    >>> mcmc.write_to_path("pythonidae.mcmc.newick", "newick")

Converting a single tree from Newick schema to NEXUS::

    >>> import dendropy
    >>> mle = dendropy.Tree.get_from_path("pythonidae.mle.newick", "newick")
    >>> mle.write_to_path("pythonidae.mle.nex", "nexus")

Collecting data from multiple sources and writing to a NEXUS-formatted file::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus", taxon_set=ds.taxon_sets[0])
    >>> ds.write_to_path("pythonidae_combined.nex", "nexus")

Note how, after the first data source has been loaded, the resulting |TaxonSet| (i.e., the first one) is passed to the subsequent :meth:`read_from_path()` statements, to ensure that the same taxa are referenced as objects corresponding to the additional data sources are created. Otherwise, as each data source is read, a new |TaxonSet| will be created, and this will result in multiple |TaxonSet| objects in the |DataSet|, with the data from each data source associated with their own, distinct |TaxonSet|.

A better way to do this is to use the "attached taxon set" mode |DataSet| object::

    >>> import dendropy
    >>> ds = dendropy.DataSet(attached_taxon_set=True)
    >>> ds.read_from_path("pythonidae_cytb.fasta", "dnafasta")
    >>> ds.read_from_path("pythonidae_aa.nex", "nexus")
    >>> ds.read_from_path("pythonidae_morphological.nex", "nexus")
    >>> ds.read_from_path("pythonidae.mle.tre", "nexus")
    >>> ds.write_to_path("pythonidae_combined.nex", "nexus")




