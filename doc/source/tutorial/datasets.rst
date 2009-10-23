********
Datasets
********

The :class:`dendropy.datasets.Dataset` class is the primary data object in |DendroPy|. It has three main attributes:

    * :attr:`taxa_blocks`
    * :attr:`trees_blocks`
    * :attr:`char_blocks`
    
which are lists of :attr:`dendropy.taxa.TaxaBlock`, :attr:`dendropy.trees.Tree` and :attr:`dendropy.characters.CharactersBlock` objects respectively.

.. toctree::
    :maxdepth: 2

    reading_and_writing.rst    
    accessing_data.rst
    