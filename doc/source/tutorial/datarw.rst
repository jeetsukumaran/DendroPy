**********************************
Creating, Reading and Writing Data
**********************************

Trees
=====

Creating and Reading Trees
--------------------------

Creating a New Tree
^^^^^^^^^^^^^^^^^^^

Phylogenetic trees in DendroPy are represented by the :class:`~dendropy.dataobject.Tree` class.
As with all other DendroPy phylogenetic data objects, this name is imported into the :mod:`dendropy` namespace, so, to create an empty :class:`~dendropy.dataobject.Tree`, you would::

    >>> import dendropy
    >>> tree = dendropy.Tree()

Or::

    >>> from dendropy import Tree
    >>> tree = Tree()

You can get some information regarding this new instance of :class:`~dendropy.dataobject.Tree` by calling the method :meth:`~dendropy.dataobject.Tree.description()`::

    >>> print(tree.description())
    Tree object at 0x79ad0 (Tree498384): ()

As you can see, it is at present an empty object. DendroPy provides a rich suite of methods to build up a tree programmatically, i.e., by creating child nodes and adding them to the tree in the proper places, setting branch length values, assigning taxa, etc.
But more than likely you would probably prefer to create a tree from a data source such as a NEXUS or NEWICK file.

If you know the data source of the tree at the time of instantiation of a new :class:`~dendropy.dataobject.Tree` object, you can pass a file-like object and a format specification string ("nexus", "newick", "nexml", etc.) to the constructor using the `stream` and `format` keyword arguments, respectively::

    >>> from dendropy import Tree
    >>> tree = Tree(stream=open('pythonidae.mcmc.nex'), format='nexus')
    >>> print(tree.description())
    Tree object at 0x553b50 (Tree6688336: 'rep.1'): (((((((('Python brongersmai':0.1,'Morelia carinata':0.1):0.1,'Morelia oenpelliensis':0.1):0.1,'Bothrochilus boa':0.1):0.1,('Antaresia perthensis':0.1,('Antaresia stimsoni':0.1,'Antaresia maculosa':0.1):0.1):0.103981):0.046169,('Python timoriensis':0.1,'Morelia bredli':0.1):0.161794):0.1,'Liasis fuscus':0.1):0.1,'Morelia boeleni':0.1):0.1,'Morelia viridis':0.1,'Aspidites ramsayi':0.1)

If the file has multiple trees (as does the :download:`the example file </examples/pythonidae.mcmc.nex>` above), then you can use the `from_index` argument to specify a *0-based index* of the particular tree that you want (i.e., the first tree in the file has an index of 0, the second has an index of 1, etc.). For example, to get the first post-burn-in tree from the :download:`the example file </examples/pythonidae.mcmc.nex>`, you could:

    >>> from dendropy import Tree
    >>> tree = Tree(stream=open('pythonidae.mcmc.nex'), format='nexus', from_index=200)
    >>> print(tree.description())
    Tree object at 0x553af0 (Tree6428592: 'rep.200000'): ('Morelia boeleni':0.106115,((('Bothrochilus boa':0.092919,'Python timoriensis':0.180712):0.020687,(((('Antaresia perthensis':0.167512,'Antaresia stimsoni':0.059787):0.033053,'Antaresia maculosa':0.146173):0.016954,(('Morelia carinata':0.100305,'Morelia bredli':0.114501):0.015794,'Morelia viridis':0.130131):0.004453):0.033047,'Liasis fuscus':0.166956):0.026128):0.004973,('Morelia oenpelliensis':0.084937,'Python brongersmai':0.245248):0.017803):0.030474,'Aspidites ramsayi':0.121686)

To build a new :class:`~dendropy.dataobject.Tree` object from a string::

    >>> from StringIO import StringIO
    >>> from dendropy import Tree
    >>> s = StringIO('((A,B),(C,D))')
    >>> tree = Tree(stream=s, format='newick')
    >>> print(tree.description())
    Tree object at 0x660e50 (Tree6339376): ((A,B),(C,D))

For convenience, you can use one of the special factory class methods of :class:`~dendropy.dataobject.Tree` to wrap up the construction. All these methods take two arguments --- an object specifying or representing the source data, and a format specification string ("nexus", "newick", "nexml", etc) specifying the format of the data --- and return a :class:`~dendropy.dataobject.Tree` object corresponding to the specified data source. They each take different types of source objects, though.
:meth:`~dendropy.dataobject.Tree.get_from_string()`, takes a string containing the data to be read as its source, :meth:`~dendropy.dataobject.Tree.get_from_path()` takes a string specifying the path to a file on the filesystem as its source, while :meth:`~dendropy.dataobject.Tree.get_from_file()` takes a file or file-like object as its source. For example::

    >>> from dendropy import Tree
    >>> tree1 = Tree.get_from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree.get_from_path('pythonidae.mcmc.nex', 'nexus')
    >>> tree3 = Tree.get_from_file(open('pythonidae.mcmc.nex', 'ru'), 'nexus')
    >>> tree4 = Tree.get_from_path('pythonidae.mcmc.nex', 'nexus', from_index=201)

(Re)defining an Existing Tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you already have an existing :class:`~dendropy.dataobject.Tree` object, and you want to redefine it or repopulate with new data, you would call one of its "read" methods:

    - :meth:`~dendropy.dataobject.Tree.read_from_file()`
    - :meth:`~dendropy.dataobject.Tree.read_from_path()`
    - :meth:`~dendropy.dataobject.Tree.read_from_string()`

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
    >>> tree.read_from_file(f, 'nexus')

Cloning an Existing Tree
^^^^^^^^^^^^^^^^^^^^^^^^

Finally, it is also possible to clone a :class:`~dendropy.dataobject.Tree` by passing it as an argument to the constructor

    >>> from dendropy import Tree
    >>> tree1 = Tree.get_from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree(tree1)
    >>> print(tree1.description())
    Tree object at 0x5e8550 (Tree6339408): ((A,B),(C,D))
    >>> print(tree2.description())
    Tree object at 0x60bbd0 (Tree6339824): ((A,B),(C,D))

This creates a *deep-copy* of the `tree1` and assigns it to `tree2`. Note that while the tree structural elements (i.e., the :class:`~dendropy.dataobject.Node` and :class:`~dendropy.dataobject.Edge` objects that make up a :class:`~dendropy.dataobject.Tree` object) are copied fully, the :class:`~dendropy.dataobject.TaxonSet` and :class:`~dendropy.dataobject.Taxon` objects are not.
This is evident when viewing more in-depth descriptions of the two :class:`~dendropy.dataobject.Tree` objects::

    >>> print(tree1.description(3))
    Tree object at 0x5e8550 (Tree6339408)
        [Taxon Set]
            TaxonSet object at 0x5508d0 (TaxonSet5572816): 4 Taxa
                [1/4] Taxon object at 0x60bf90 (Taxon6340496): 'A'
                [2/4] Taxon object at 0x60be30 (Taxon6340144): 'B'
                [3/4] Taxon object at 0x60bfb0 (Taxon6340528): 'C'
                [4/4] Taxon object at 0x60bcd0 (Taxon6339792): 'D'
        [Tree]
            ((A,B),(C,D))
    >>> print(tree2.description(3))
    Tree object at 0x60bbd0 (Tree6339824)
        [Taxon Set]
            TaxonSet object at 0x5508d0 (TaxonSet5572816): 4 Taxa
                [1/4] Taxon object at 0x60bf90 (Taxon6340496): 'A'
                [2/4] Taxon object at 0x60be30 (Taxon6340144): 'B'
                [3/4] Taxon object at 0x60bfb0 (Taxon6340528): 'C'
                [4/4] Taxon object at 0x60bcd0 (Taxon6339792): 'D'
        [Tree]
            ((A,B),(C,D))




This is behavior is by design, based on the logic that while you wanted a copy of the tree, you still mean them to refer to the same taxa.
So, for example, if you were to prune or move an edge, change the edge lengths, etc. on `tree2`, or even reassign a particular :class:`~dendropy.dataobject.Taxon` object to a different node, it would not in any way affect `tree1`.
But if you were to assign a different label to a :class:`~dendropy.dataobject.Taxon` object on `tree2`, this *would* affect `tree11`.
