*******************************
Phylogenetic Character Analyses
*******************************

Phylogenetic Independent Contrasts (PIC)
========================================

Basic Analysis
--------------

A phylogenetic independent contrasts analysis (Felsenstein 1985; Garland et al. 2005) can be carried out using the :class:`~dendropy.continuous.PhylogeneticIndependentConstrasts` class.
This requires you to have a |Tree| and a |ContinuousCharacterMatrix| which reference the same |TaxonSet|.
Thus, if your data is in the same file::

    >>> import dendropy
    >>> dataset = dendropy.DataSet.get_from_path("primates.cont.nex", "nexus")
    >>> tree = dataset.tree_list[0][0]
    >>> chars = dataset.char_matrices[0]

While if you have the tree and characters in a different file::

    >>> import dendropy
    >>> taxa = dendropy.TaxonSet()
    >>> tree = dendropy.Tree.get_from_path("primates.tre", "newick", taxon_set=taxa)
    >>> chars = dendropy.ContinuousCharacterMatrix.get_from_path("primates.cc.nex", "nexus", taxon_set=taxa)

In either case, we have a |Tree| object, ``tree`` and a |ContinuousCharacterMatrix| object, ``chars``, that both reference the same |TaxonSet|.

Once the data is loaded, we create the :class:`~dendropy.continuous.PhylogeneticIndependentConstrasts` object::

    >>> from dendropy import continuous
    >>> pic = dendropy.continuous.PhylogeneticIndependentContrasts(tree=tree, char_matrix=chars)

At this point, the data is ready for analysis. Typically, we want to map the contrasts onto a tree. The :meth:`~dendropy.continuous.PhylogeneticIndependentConstrasts.contrasts_tree` method takes a single mandatory argument, the 0-based index of the character (or column) to be analyzed, and returns a |Tree| object that is a clone of the original input |Tree|, but with the following attributes added to each |Node|:

        - `pic_state_value`
        - `pic_state_variance`
        - `pic_contrast_raw`
        - `pic_contrast_variance`
        - `pic_contrast_standardized`
        - `pic_edge_length_error`
        - `pic_corrected_edge_length`

In addition to the 0-based index first argument, ``character_index``, the :meth:`~dendropy.continuous.PhylogeneticIndependentConstrasts.contrasts_tree` method takes the following optional arguments:

    ``annotate_pic_statistics``
        If |True| then the PIC statistics attributes will be *annotated* (i.e., serialized or persisted when the tree is written out or saved. Defaults to |False|.
    ``state_values_as_node_labels``
        If |True| then the :class:`~dendropy.dataobject.Tree.Node.label` attribute of each |Node| object will be set to the value of the character.
    ``corrected_edge_lengths``
        If |True| then the |Tree| returned will have its edge lengths adjusted to the corrected edge lengths as yielded by the PIC analysis.

Results as a Table
------------------

So the following retrieves the constrasts tree for the first character (index=0), and prints a table of the various statistics::

    >>> ctree1 = pic.contrasts_tree(character_index=0,
    ...     annotate_pic_statistics=True,
    ...     state_values_as_node_labels=False,
    ...     corrected_edge_lengths=False)
    >>> for nd in ctree1.postorder_internal_node_iter():
    ...     row = [nd.pic_state_value,
    ...             nd.pic_state_variance,
    ...             nd.pic_contrast_raw,
    ...             nd.pic_edge_length_error]
    ...     row_str = [(("%10.8f") % i) for i in row]
    ...     row_str = "    ".join(row_str)
    ...     label = nd.label.ljust(6)
    ...     print "%s %s" % (label, row_str)
    HP     3.85263000    0.38500000    0.48342000    0.10500000
    HPM    3.20037840    0.34560000    1.48239000    0.21560000
    HPMA   2.78082358    0.60190555    1.17222840    0.22190555
    Root   1.18372461    0.37574347    4.25050358    0.37574347

Results as a Newick String with State Values as Node Labels
-----------------------------------------------------------

Alternatively, you might want to visualize the results as a tree showing the numeric values of the states. The following produces this for each character in the matrix by first requesting that :meth:`~dendropy.continuous.PhylogeneticIndependentConstrasts.contrasts_tree` replace existing node labels with the state values for that node, and then, when writing out in Newick format, suppressing taxon labels and printing node labels in their place:

.. literalinclude:: /examples/pic1.py

This results in::

    [&R] ((((4.09434:0.21,3.61092:0.21)3.85263:0.28,2.37024:0.49)3.2003784:0.13,2.02815:0.62)2.78082357912:0.38,'-1.46968':1.0)1.1837246134:0.0;

    [&R] ((((4.74493:0.21,3.3322:0.21)4.038565:0.28,3.3673:0.49)3.7432084:0.13,2.89037:0.62)3.43796714996:0.38,2.30259:1.0)3.01135659943:0.0;

Results as a NEXUS Document with Analysis Statistics as Node Metadata
---------------------------------------------------------------------

However, probably the best way to visualize the results would be as a tree marked up with metadata that can be viewed in |FigTree|_ (by checking "Node Labels" and selecting the appropriate statistics from the drop-down menu). This is, in fact, even easier to do than the above, as it will result from the default options. The following illustrates this. It collects the metadata-annotated contrast analysis trees produced by :meth:`~dendropy.continuous.PhylogeneticIndependentConstrasts.contrasts_tree` in a |TreeList| object, and then prints the |TreeList| as NEXUS-formatted string. The default options to :meth:`~dendropy.continuous.PhylogeneticIndependentConstrasts.contrasts_tree` result in annotated attributes, while the default options to the writing method result in the annotations being written out as comment metadata.

.. literalinclude:: /examples/pic_annotated.py

Thus, we get::

    #NEXUS


    BEGIN TAXA;
        DIMENSIONS NTAX=5;
        TAXLABELS
            Homo
            Pongo
            Macaca
            Ateles
            Galago
    ;
    END;

    BEGIN TREES;
        TREE PIC_1 = [&R] ((((Homo:0.21[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.21,pic_state_value=4.09434,pic_contrast_standardized=None,pic_contrast_raw=None],Pongo:0.21[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.21,pic_state_value=3.61092,pic_contrast_standardized=None,pic_contrast_raw=None])HP:0.28[&pic_contrast_variance=0.42,pic_edge_length_error=0.105,pic_state_variance=0.385,pic_corrected_edge_length=0.385,pic_state_value=3.85263,pic_contrast_standardized=0.745933254387,pic_contrast_raw=0.48342],Macaca:0.49[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.49,pic_state_value=2.37024,pic_contrast_standardized=None,pic_contrast_raw=None])HPM:0.13[&pic_contrast_variance=0.875,pic_edge_length_error=0.2156,pic_state_variance=0.3456,pic_corrected_edge_length=0.3456,pic_state_value=3.2003784,pic_contrast_standardized=1.58474156959,pic_contrast_raw=1.48239],Ateles:0.62[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.62,pic_state_value=2.02815,pic_contrast_standardized=None,pic_contrast_raw=None])HPMA:0.38[&pic_contrast_variance=0.9656,pic_edge_length_error=0.221905550953,pic_state_variance=0.601905550953,pic_corrected_edge_length=0.601905550953,pic_state_value=2.78082357912,pic_contrast_standardized=1.19292629182,pic_contrast_raw=1.1722284],Galago:1.0[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=1.0,pic_state_value=-1.46968,pic_contrast_standardized=None,pic_contrast_raw=None])Root:0.0[&pic_contrast_variance=1.60190555095,pic_edge_length_error=0.37574347039,pic_state_variance=0.37574347039,pic_corrected_edge_length=0.37574347039,pic_state_value=1.1837246134,pic_contrast_standardized=3.35831889583,pic_contrast_raw=4.25050357912];

        TREE PIC_2 = [&R] ((((Homo:0.21[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.21,pic_state_value=4.74493,pic_contrast_standardized=None,pic_contrast_raw=None],Pongo:0.21[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.21,pic_state_value=3.3322,pic_contrast_standardized=None,pic_contrast_raw=None])HP:0.28[&pic_contrast_variance=0.42,pic_edge_length_error=0.105,pic_state_variance=0.385,pic_corrected_edge_length=0.385,pic_state_value=4.038565,pic_contrast_standardized=2.17988971592,pic_contrast_raw=1.41273],Macaca:0.49[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.49,pic_state_value=3.3673,pic_contrast_standardized=None,pic_contrast_raw=None])HPM:0.13[&pic_contrast_variance=0.875,pic_edge_length_error=0.2156,pic_state_variance=0.3456,pic_corrected_edge_length=0.3456,pic_state_value=3.7432084,pic_contrast_standardized=0.717612470209,pic_contrast_raw=0.671265],Ateles:0.62[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=0.62,pic_state_value=2.89037,pic_contrast_standardized=None,pic_contrast_raw=None])HPMA:0.38[&pic_contrast_variance=0.9656,pic_edge_length_error=0.221905550953,pic_state_variance=0.601905550953,pic_corrected_edge_length=0.601905550953,pic_state_value=3.43796714996,pic_contrast_standardized=0.867896862104,pic_contrast_raw=0.8528384],Galago:1.0[&pic_contrast_variance=None,pic_edge_length_error=0.0,pic_state_variance=None,pic_corrected_edge_length=1.0,pic_state_value=2.30259,pic_contrast_standardized=None,pic_contrast_raw=None])Root:0.0[&pic_contrast_variance=1.60190555095,pic_edge_length_error=0.37574347039,pic_state_variance=0.37574347039,pic_corrected_edge_length=0.37574347039,pic_state_value=3.01135659943,pic_contrast_standardized=0.897060422518,pic_contrast_raw=1.13537714996];

    END;


Multifurcating Trees and Polytomies
-----------------------------------

By default, the :class:`~dendropy.continuous.PhylogeneticIndependentConstrasts` class only handles fully-bifurcating trees, and throws an exception if the input tree has polytomies.
You can change this behavior by specifying one of the following strings to the "``polytomy_strategy``" argument of the class constructor:

    "``ignore``"
        Polytomies will handled without complaint::

            >>> pic = dendropy.continuous.PhylogeneticIndependentContrasts(tree=tree,
            ...        char_matrix=chars,
            ...        polytomy_strategy='ignore')

        Note that in this case the raw contrast and the raw contrast variance calculated for nodes that have more than two children will be invalid. The reconstructed state values should be still valid, though.

    "``resolve``"
        Polytomies will be arbitrarily resolved with 0-length branches::

            >>> pic = dendropy.continuous.PhylogeneticIndependentContrasts(tree=tree,
            ...        char_matrix=chars,
            ...        polytomy_strategy='resolve')

        In this case this validity of the analysis for nodes with (originally) more than two children is dubious, as the resulting contrasts are non-independent.



