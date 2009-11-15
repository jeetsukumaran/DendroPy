**********************************
Creating, Reading and Writing Data
**********************************

Trees
=====

Creating and Reading Trees
--------------------------

Phylogenetic trees in DendroPy are represented by the :class:`~dendropy.dataobject.Tree` class.
As with all other DendroPy phylogenetic data objects, this name is imported into the :mod:`dendropy` namespace, so, to create an empty :class:`~dendropy.dataobject.Tree`, you would::

    >>> import dendropy
    >>> tree = dendropy.Tree()

Or::

    >>> from dendropy import Tree
    >>> tree = Tree()

You can get some information regarding this new instance of :class:`~dendropy.dataobject.Tree` by calling the method :meth:`~dendropy.dataobject.Tree.describe()`::

    >>> tree.describe()

As you can see, it is at present an empty object. DendroPy provides a rich suite of methods to build up a tree programmatically, i.e., by creating child nodes and adding them to the tree in the proper places, setting branch length values, assigning taxa, etc.
But more than likely you would probably prefer to create a tree from a data source such as a NEXUS or NEWICK file.

If you know the data source of the tree at the time of instantiation of a new :class:`~dendropy.dataobject.Tree` object, you can pass a file-like object and a format specification string ("nexus", "newick", "nexml", etc.) to the constructor using the `stream` and `format` keyword arguments, respectively::

    >>> from dendropy import Tree
    >>> tree = Tree(stream=open('mle.tre'), format='newick')

Or::

    >>> from StringIO import StringIO
    >>> from dendropy import Tree
    >>> s = StringIO('((A,B),(C,D))')
    >>> tree = Tree(stream=s, format='newick')

For convenience, you can use one of the special factory class methods of :class:`~dendropy.dataobject.Tree` to wrap up the construction. All these methods take two arguments --- an object specifying or representing the source data, and a format specification string ("nexus", "newick", "nexml", etc) specifying the format of the data --- and return a :class:`~dendropy.dataobject.Tree` object corresponding to the specified data source. They each take different types of source objects, though.
:meth:`~dendropy.dataobject.Tree.from_string()`, takes a string containing the data to be read as its source, :meth:`~dendropy.dataobject.Tree.from_path()` takes a string specifying the path to a file on the filesystem as its source, while :meth:`~dendropy.dataobject.Tree.from_file()` takes a file or file-like object as its source. For example::

    >>> from dendropy import Tree
    >>> tree1 = Tree.from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree.from_path('mle.tre', 'newick')
    >>> tree3 = Tree.from_file(open('mle.tre', 'ru'), 'newick')

If you already have a :class:`~dendropy.dataobject.Tree`, and you want to redefine it or repopulate with new data, you would call one of its "read" methods:

    - :meth:`~dendropy.dataobject.Tree.read_file()`
    - :meth:`~dendropy.dataobject.Tree.read_path()`
    - :meth:`~dendropy.dataobject.Tree.read_string()`

For example::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> tree.read_string('((A,B),(C,D))', 'newick')

Or reading from a file path::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> tree.read_path('mle.tre', 'newick')

Or a file object::

    >>> from dendropy import Tree
    >>> tree = Tree()
    >>> f = open('mle.tre', 'rU')
    >>> tree.read_file(f, 'newick')

Finally, it is also possible to clone a :class:`~dendropy.dataobject.Tree` by passing it as an argument to the constructor

    >>> from dendropy import Tree
    >>> tree1 = Tree.from_string('((A,B),(C,D))', 'newick')
    >>> tree2 = Tree(tree1)
