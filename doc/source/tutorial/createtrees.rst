**************************
Creating and Reading Trees
**************************

Creating a New Tree
===================

Phylogenetic trees in DendroPy are represented by the :class:`~dendropy.dataobject.tree.Tree` class.
As with all other DendroPy phylogenetic data objects, this name is imported into the :mod:`dendropy` namespace, so, to create an empty :class:`~dendropy.dataobject.tree.Tree`, you would::

    >>> import dendropy
    >>> tree = dendropy.Tree()

Or::

    >>> from dendropy import Tree
    >>> tree = Tree()

You can get some information regarding this new instance of :class:`~dendropy.dataobject.tree.Tree` by calling the method :meth:`~dendropy.dataobject.tree.Tree.description()`::

    >>> print(tree.description())
    Tree object at 0x79ad0 (Tree498384): ()

As you can see, it is at present an empty object. DendroPy provides a rich suite of methods to build up a tree programmatically, i.e., by creating child nodes and adding them to the tree in the proper places, setting branch length values, assigning taxa, etc.
But more than likely you would probably prefer to create a tree from a data source such as a NEXUS or NEWICK file.

If you know the data source of the tree at the time of instantiation of a new :class:`~dendropy.dataobject.tree.Tree` object, you can pass a file-like object and a format specification string ("nexus", "newick", "nexml", etc.) to the constructor using the `stream` and `format` keyword arguments, respectively::

    >>> from dendropy import Tree
    >>> tree = Tree(stream=open('pythonidae.mcmc.nex'), format='nexus')
    >>> print(tree.description())
    Tree object at 0x553b50 (Tree6688336: 'rep.1'): (((((((('Python brongersmai':0.1,'Morelia carinata':0.1):0.1,'Morelia oenpelliensis':0.1):0.1,'Bothrochilus boa':0.1):0.1,('Antaresia perthensis':0.1,('Antaresia stimsoni':0.1,'Antaresia maculosa':0.1):0.1):0.103981):0.046169,('Python timoriensis':0.1,'Morelia bredli':0.1):0.161794):0.1,'Liasis fuscus':0.1):0.1,'Morelia boeleni':0.1):0.1,'Morelia viridis':0.1,'Aspidites ramsayi':0.1)

If the file has multiple trees (as does the :download:`the example file </examples/pythonidae.mcmc.nex>` above), then you can use the `from_index` argument to specify a *0-based index* of the particular tree that you want (i.e., the first tree in the file has an index of 0, the second has an index of 1, etc.). For example, to get the first post-burn-in tree from the :download:`the example file </examples/pythonidae.mcmc.nex>`, you could:

    >>> from dendropy import Tree
    >>> tree = Tree(stream=open('pythonidae.mcmc.nex'), format='nexus', from_index=200)
    >>> print(tree.description())
    Tree object at 0x553af0 (Tree6428592: 'rep.200000'): ('Morelia boeleni':0.106115,((('Bothrochilus boa':0.092919,'Python timoriensis':0.180712):0.020687,(((('Antaresia perthensis':0.167512,'Antaresia stimsoni':0.059787):0.033053,'Antaresia maculosa':0.146173):0.016954,(('Morelia carinata':0.100305,'Morelia bredli':0.114501):0.015794,'Morelia viridis':0.130131):0.004453):0.033047,'Liasis fuscus':0.166956):0.026128):0.004973,('Morelia oenpelliensis':0.084937,'Python brongersmai':0.245248):0.017803):0.030474,'Aspidites ramsayi':0.121686)

To build a new :class:`~dendropy.dataobject.tree.Tree` object from a string::

    >>> from StringIO import StringIO
    >>> from dendropy import Tree
    >>> s = StringIO('((A,B),(C,D))')
    >>> tree = Tree(stream=s, format='newick')
    >>> print(tree.description())
    Tree object at 0x660e50 (Tree6339376): ((A,B),(C,D))

For convenience, you can use one of the special factory class methods of :class:`~dendropy.dataobject.tree.Tree` to wrap up the construction. All these methods take two arguments --- an object specifying or representing the source data, and a format specification string ("nexus", "newick", "nexml", etc) specifying the format of the data --- and return a :class:`~dendropy.dataobject.tree.Tree` object corresponding to the specified data source. They each take different types of source objects, though.
:meth:`~dendropy.dataobject.tree.Tree.get_from_string()`, takes a string containing the data to be read as its source, :meth:`~dendropy.dataobject.tree.Tree.get_from_path()` takes a string specifying the path to a file on the filesystem as its source, while :meth:`~dendropy.dataobject.tree.Tree.get_from_stream()` takes a file or file-like object as its source. For example::

    >>> from dendropy import Tree
    >>> tree1 = Tree.get_from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree.get_from_path('pythonidae.mcmc.nex', 'nexus')
    >>> tree3 = Tree.get_from_stream(open('pythonidae.mcmc.nex', 'ru'), 'nexus')
    >>> tree4 = Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=201)

Reading into an Existing Tree
=============================

If you already have an existing :class:`~dendropy.dataobject.tree.Tree` object, and you want to redefine it or repopulate with new data, you would call one of its "read" methods:

    - :meth:`~dendropy.dataobject.tree.Tree.read_from_stream()`
    - :meth:`~dendropy.dataobject.tree.Tree.read_from_path()`
    - :meth:`~dendropy.dataobject.tree.Tree.read_from_string()`

For example::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> tree.read_from_string('((A,B),(C,D))', 'newick')

Or reading from a file path::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> tree.read_from_path('pythonidae.mcmc.nex', 'nexus')

Or a file object::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> f = open('pythonidae.mcmc.nex', 'rU')
    >>> tree.read_from_stream(f, 'nexus')

Cloning an Existing Tree
========================

Finally, it is also possible to clone a :class:`~dendropy.dataobject.tree.Tree` by passing it as an argument to the constructor

    >>> from dendropy import Tree
    >>> tree1 = Tree.get_from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree(tree1)
    >>> print(tree1.description())
    Tree object at 0x5e8550 (Tree6339408): ((A,B),(C,D))
    >>> print(tree2.description())
    Tree object at 0x60bbd0 (Tree6339824): ((A,B),(C,D))

This creates a *deep-copy* of the `tree1` and assigns it to `tree2`. Note that while the tree structural elements (i.e., the :class:`~dendropy.dataobject.tree.Node` and :class:`~dendropy.dataobject.tree.Edge` objects that make up a :class:`~dendropy.dataobject.tree.Tree` object) are copied fully, the :class:`~dendropy.dataobject.taxon.TaxonSet` and :class:`~dendropy.dataobject.taxon.Taxon` objects are not.
This is evident when viewing more in-depth descriptions of the two :class:`~dendropy.dataobject.tree.Tree` objects::

    >>> print(tree1.description(3))
    Tree object at 0x79bd0 (Tree6679376): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonSet object at 0x7c7b0 (TaxonSet509872): 4 Taxa
                [1/4] Taxon object at 0x65ec90 (Taxon6679696): 'A'
                [2/4] Taxon object at 0x65ecf0 (Taxon6679792): 'B'
                [3/4] Taxon object at 0x65ed90 (Taxon6679952): 'C'
                [4/4] Taxon object at 0x65edf0 (Taxon6680048): 'D'
        [Tree]
            ((A,B),(C,D))
        [Nodes]
            [1/7] Node object at 0x65ebd0 (Node6679504)
            [2/7] Node object at 0x65ec10 (Node6679568)
            [3/7] Node object at 0x65ec50 (Node6679632)
            [4/7] Node object at 0x65ecb0 (Node6679728)
            [5/7] Node object at 0x65ed10 (Node6679824)
            [6/7] Node object at 0x65ed50 (Node6679888)
            [7/7] Node object at 0x65edb0 (Node6679984)
        [Edges]
            [1/7] Edge object at 0x65ebf0 (Edge6679536, Length=None)
            [2/7] Edge object at 0x65ec30 (Edge6679600, Length=None)
            [3/7] Edge object at 0x65ec70 (Edge6679664, Length=None)
            [4/7] Edge object at 0x65ecd0 (Edge6679760, Length=None)
            [5/7] Edge object at 0x65ed30 (Edge6679856, Length=None)
            [6/7] Edge object at 0x65ed70 (Edge6679920, Length=None)
            [7/7] Edge object at 0x65edd0 (Edge6680016, Length=None)
    >>> print(tree2.description(3))
    Tree object at 0x65eb50 (Tree6680080): 7 Nodes, 7 Edges
        [Taxon Set]
            TaxonSet object at 0x7c7b0 (TaxonSet509872): 4 Taxa
                [1/4] Taxon object at 0x65ec90 (Taxon6679696): 'A'
                [2/4] Taxon object at 0x65ecf0 (Taxon6679792): 'B'
                [3/4] Taxon object at 0x65ed90 (Taxon6679952): 'C'
                [4/4] Taxon object at 0x65edf0 (Taxon6680048): 'D'
        [Tree]
            ((A,B),(C,D))
        [Nodes]
            [1/7] Node object at 0x65ee90 (Node6680208)
            [2/7] Node object at 0x65eeb0 (Node6680240)
            [3/7] Node object at 0x65ef10 (Node6680336)
            [4/7] Node object at 0x65ef90 (Node6680464)
            [5/7] Node object at 0x65ef30 (Node6680368)
            [6/7] Node object at 0x65efd0 (Node6680528)
            [7/7] Node object at 0x64b070 (Node6598768)
        [Edges]
            [1/7] Edge object at 0x65eef0 (Edge6680304, Length=None)
            [2/7] Edge object at 0x65ef50 (Edge6680400, Length=None)
            [3/7] Edge object at 0x65efb0 (Edge6680496, Length=None)
            [4/7] Edge object at 0x65eff0 (Edge6680560, Length=None)
            [5/7] Edge object at 0x64b030 (Edge6598704, Length=None)
            [6/7] Edge object at 0x64b090 (Edge6598800, Length=None)
            [7/7] Edge object at 0x64b0d0 (Edge6598864, Length=None)

As you can see, the :class:`~dendropy.dataobject.tree.Node` and :class:`~dendropy.dataobject.tree.Edge` objects are distinct between the trees, but the associated taxa and taxon references are the same.
This is based on the logic that while you want an independent copy of the tree, you still dealing with the same taxa.
So, for example, if you were to prune or move an edge, change the edge lengths, etc. on `tree2`, or even reassign a particular :class:`~dendropy.dataobject.taxon.Taxon` object to a different node, it would not in any way affect `tree1`.
But if you were to assign a different label to a :class:`~dendropy.dataobject.taxon.Taxon` object on `tree2`, this *would* affect the same :class:`~dendropy.dataobject.taxon.Taxon` object on `tree11`.

Taxon Management
================
Every time an independent :class:`~dendropy.dataobject.tree.Tree` object is created, by default a new :class:`~dendropy.dataobject.taxon.TaxonSet` object is created and associated with the :class:`~dendropy.dataobject.tree.Tree`. 
This means that if two :class:`~dendropy.dataobject.tree.Tree` objects are independentally created, even if from the same data source, they will reference distinct sets of :class:`~dendropy.dataobject.taxon.Taxon` objects (though the labels might be the same).
Consider the following::

    >>> t1 = dendropy.Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=199)
    >>> t2 = dendropy.Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=200)
    
Here, two tree objects are created from the 200th and 201st trees defined in the :download:`the example file </examples/pythonidae.mcmc.nex>`. 
Even though they were sourced from the same data file, and, indeed, the same "TREE" block within the same data file, as a result of their independent default instantiation, they refer to distinct (though similar) :class:`~dendropy.dataobject.taxon.TaxonSet` and :class:`~dendropy.dataobject.taxon.Taxon` objects::

    >>> print(t1.taxon_set.description(2))
    TaxonSet object at 0x101f630 (TaxonSet16905776): 13 Taxa
        [0] Taxon object at 0x124e690 (Taxon19195536): 'Aspidites ramsayi'
        [1] Taxon object at 0x124e630 (Taxon19195440): 'Bothrochilus boa'
        [2] Taxon object at 0x124e090 (Taxon19194000): 'Liasis fuscus'
        [3] Taxon object at 0x124e670 (Taxon19195504): 'Antaresia stimsoni'
        [4] Taxon object at 0x124e330 (Taxon19194672): 'Morelia viridis'
        [5] Taxon object at 0x124e6d0 (Taxon19195600): 'Morelia bredli'
        [6] Taxon object at 0x124e6b0 (Taxon19195568): 'Antaresia perthensis'
        [7] Taxon object at 0x124e3d0 (Taxon19194832): 'Python timoriensis'
        [8] Taxon object at 0x124e6f0 (Taxon19195632): 'Antaresia maculosa'
        [9] Taxon object at 0x124e350 (Taxon19194704): 'Morelia carinata'
        [10] Taxon object at 0x124e710 (Taxon19195664): 'Python brongersmai'
        [11] Taxon object at 0x124e650 (Taxon19195472): 'Morelia boeleni'
        [12] Taxon object at 0x124e770 (Taxon19195760): 'Morelia oenpelliensis'
    >>> print(t2.taxon_set.description(2))
    TaxonSet object at 0x1243a80 (TaxonSet19151488): 13 Taxa
        [0] Taxon object at 0x129e610 (Taxon19523088): 'Aspidites ramsayi'
        [1] Taxon object at 0x129e930 (Taxon19523888): 'Bothrochilus boa'
        [2] Taxon object at 0x129e7b0 (Taxon19523504): 'Liasis fuscus'
        [3] Taxon object at 0x129e850 (Taxon19523664): 'Antaresia stimsoni'
        [4] Taxon object at 0x129e950 (Taxon19523920): 'Morelia viridis'
        [5] Taxon object at 0x129e770 (Taxon19523440): 'Morelia bredli'
        [6] Taxon object at 0x129e4b0 (Taxon19522736): 'Antaresia perthensis'
        [7] Taxon object at 0x129ec70 (Taxon19524720): 'Python timoriensis'
        [8] Taxon object at 0x129eb90 (Taxon19524496): 'Antaresia maculosa'
        [9] Taxon object at 0x129e7d0 (Taxon19523536): 'Morelia carinata'
        [10] Taxon object at 0x129e550 (Taxon19522896): 'Python brongersmai'
        [11] Taxon object at 0x129e730 (Taxon19523376): 'Morelia boeleni'
        [12] Taxon object at 0x129e830 (Taxon19523632): 'Morelia oenpelliensis'
        
This would mean that any comparisons between the two trees would be invalid::

    >>> from dendropy import treecalc
    >>> treecalc.symmetric_difference(t1, t2)
    ------------------------------------------------------------
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    
      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 280, in symmetric_difference
        t = false_positives_and_negatives(tree1, tree2)
    
      File "/Users/jeet/Documents/Projects/dendropy/dendropy/treecalc.py", line 293, in false_positives_and_negatives
        % (hex(id(reference_tree.taxon_set)), hex(id(test_tree.taxon_set))))
    
    TypeError: Trees have different TaxonSet objects: 0x101f630 vs. 0x1243a80
    
At this stage, short of r    
    
The correct way to instantiate two :class:`~dendropy.dataobject.tree.Tree` objects so that they refer to the same taxa objects is to pass a :class:`~dendropy.dataobject.taxon.TaxonSet` for them to use::

    >>> t1 = dendropy.Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=199)
    >>> t2 = dendropy.Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=200, taxon_set=t1.taxon_set)
    >>> treecalc.symmetric_difference(t1, t2)
    8    

The same applies even if they are sourced from different files: specifying a :class:`~dendropy.dataobject.taxon.TaxonSet` object ensures that the different :class:`~dendropy.dataobject.tree.Tree` objects reference the same taxa::

    >>> t1 = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
    >>> t2 = dendropy.Tree.get_from_path('pythonidae.pars.newick', 'newick', taxon_set=t1.taxon_set)
    >>> treecalc.symmetric_difference(t1, t2)
    4
    
It probably would lead to more maintainable code if you were to explicitly create a :class:`~dendropy.dataobject.taxon.TaxonSet` object, and pass that to all :class:`~dendropy.dataobject.tree.Tree` objects that you create (assuming, of course, that they all do indeed refer to the same taxa)::

    >>> taxa = dendropy.TaxonSet()
    >>> t1 = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus', taxon_set=taxa) 
    >>> t2 = dendropy.Tree.get_from_path('pythonidae.pars.newick', 'newick', taxon_set=taxa) 
    >>> t3 = dendropy.Tree.get_from_path('pythonidae.mcmc.nexus', 'nexus', from_index=199, taxon_set=taxa)
    
There is no doubt that, while this approach is acceptable for a small number of specific trees, it can be tedious and error-prone as things scale up.
The preferred way of dealing with multiple trees referencing the same taxa is to use a :class:`~dendropy.dataobject.tree.TreeList` object to collect and manage them, as this automatically enforces the homogeneity of :class:`~dendropy.dataobject.taxon.TaxonSet` references amongst all its members.
This is covered in the :doc:`next section <createtreelists>`
