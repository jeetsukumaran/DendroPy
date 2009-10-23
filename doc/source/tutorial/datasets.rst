********
Datasets
********

The :class:`dendropy.datasets.Dataset` class is the primary data management object in |DendroPy|. It has three main attributes:

    * :attr:`taxa_blocks`
    * :attr:`trees_blocks`
    * :attr:`char_blocks`
    
which are lists of :class:`dendropy.taxa.TaxaBlock`, :class:`dendropy.trees.Tree` and :class:`dendropy.characters.CharactersBlock` objects respectively.

You will use objects of the :class:`Dataset` to read, access, manage and write phylogenetic data such as taxa, trees and characters.

.. toctree::
    :maxdepth: 2

    reading_and_writing.rst    
    accessing_data.rst
    