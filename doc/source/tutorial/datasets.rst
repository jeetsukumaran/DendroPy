********
Datasets
********

The :class:`dendropy.datasets.Dataset` class is the primary data manager in |DendroPy|. It has three main attributes:

    * :attr:`taxa_blocks`
    * :attr:`trees_blocks`
    * :attr:`char_blocks`
    
:attr:`taxa_blocks` and :attr:`trees_blocks` are specialized Python lists of :class:`TaxaBlock` and :class:`TreesBlock` objects, respectively, while :attr:`char_blocks` is a :class:`dendropy.characters.CharactersBlock` object.

.. toctree::
    :maxdepth: 2

    reading_and_writing.rst    
    accessing_data.rst
    