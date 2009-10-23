***************************************
:mod:`datasets` -- Primary Data Manager
***************************************

.. module:: datasets

.. toctree::
    :maxdepth: 2

.. autoclass:: dendropy.datasets.Dataset
    :members:
    
    The :class:`dendropy.datasets.Dataset` is composed of three basic list objects:
    
        * .. attribute:: taxa_blocks
        
                A list of :class:`TaxaBlock` objects.
        
        * .. attribute:: trees_blocks
        
                A list of :class:`TreesBlock` objects.    
            
        * .. attribute:: chars_blocks
        
                A list of :class:`CharactersBlock` objects.    
    
