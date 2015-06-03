*******************************
Converting Between Data Formats
*******************************

Any data in a schema that can be read by DendroPy, can be saved to files in any schema that can be written by DendroPy.
Converting data between formats is simply a matter of calling readers and writers of the appropriate type.

Converting from FASTA schema to NEXUS::

    >>> import dendropy
    >>> cytb = dendropy.DnaCharacterMatrix.get(path="pythonidae_cytb.fasta", schema="fasta")
    >>> cytb.write(path="pythonidae_cytb.nexus", schema="nexus")

Converting a collection of trees from NEXUS schema to Newick::

    >>> import dendropy
    >>> post_trees = dendropy.TreeList()
    >>> post_trees.read(
    ...         file=open("pythonidae.nex.run1.t", "r")
    ...         schema="nexus",
    ...         tree_offset=200)
    >>> post_trees.read(
    ...         path="pythonidae.nex.run2.t",
    ...         schema="nexus",
    ...         tree_offset=200)
    >>> post_trees.write(
    ...     path="pythonidae.mcmc.newick",
    ...     schema="newick")

Converting a single tree from Newick schema to NEXUS::

    >>> import dendropy
    >>> mle = dendropy.Tree.get(path="pythonidae.mle.newick", schema="newick")
    >>> mle.write(path="pythonidae.mle.nex", schema="nexus")

Collecting data from multiple sources and writing to a NEXUS-formatted file::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> tns = ds.new_taxon_namespace()
    >>> ds.attach_taxon_namespace(tns)
    >>> ds.read(
    ...     path="pythonidae_cytb.fasta",
    ...     schema="fasta",
    ...     data_type="dna")
    >>> ds.read(
    ...     path="pythonidae_aa.nex",
    ...     schema="nexus")
    >>> ds.read(
    ...     path="pythonidae_morph.nex",
    ...     schema="nexus")
    >>> ds.read(
    ...     path="pythonidae_trees.tre",
    ...     schema="newick")
    >>> ds.write(
    ...     path="pythonidae_combined.nex",
    ...     schema="nexus")

Note how we create a new |TaxonNamespace| instance using the :meth:`~dendropy.datamodel.DataSet.new_taxon_namespace` method, and then "bind" or attach it to the |DataSet| instance using the :meth:`~dendropy.datamodel.DataSet.attach_taxon_namespace()` method.
This ensures that all new data parsed by the |DataSet| instance will reference the same |TaxonNamespace| instance, i.e., all taxon labels will be mapped to the same set of |Taxon| objects.
Alternatively, we could also have explicitly passed in the |TaxonNamespace| instance to use for each reading operation::

    >>> import dendropy
    >>> ds = dendropy.DataSet()
    >>> tns = ds.new_taxon_namespace()
    >>> ds.read(
    ...     path="pythonidae_cytb.fasta",
    ...     schema="fasta",
    ...     data_type="dna",
    ...     taxon_namespace=tns)
    >>> ds.read(
    ...     path="pythonidae_aa.nex",
    ...     schema="nexus",
    ...     taxon_namespace=tns)
    >>> ds.read(
    ...     path="pythonidae_morph.nex",
    ...     schema="nexus",
    ...     taxon_namespace=tns)
    >>> ds.read(
    ...     path="pythonidae_trees.tre",
    ...     schema="newick",
    ...     taxon_namespace=tns)
    >>> ds.write(
    ...     path="pythonidae_combined.nex",
    ...     schema="nexus",
    ...     taxon_namespace=tns)
