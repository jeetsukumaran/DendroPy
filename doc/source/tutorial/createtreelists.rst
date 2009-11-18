*******************************
Creating and Reading Tree Lists
*******************************

Multiple :class:`~dendropy.dataobject.Tree` objects referencing the same taxa can be aggregated into a  :class:`~dendropy.dataobject.TreeList`.
You can build up a :class:`~dendropy.dataobject.TreeList` object by instantiating an empty one, and individually appending  :class:`~dendropy.dataobject.Tree` objects to it::

    >>> import dendropy
    >>> trees = dendropy.TreeList()
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> trees.append(t1)
    >>> trees.append(t2)
    
Note that when `t1` and `t2` were first created, they referred to distinct    

If the :class:`~dendropy.dataobject.Tree` objects exist at the time of instantiation of the :class:`~dendropy.dataobject.TreeList`, you could pass them as a list to the constructor:

    >>> import dendropy
    >>> t1 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> t2 = dendropy.Tree.get_from_string('((A,B),(C,D))', format='newick')
    >>> trees = dendropy.TreeList([t1, t2])

Alternatively, you can use one of the factory methods (:meth:`~dendropy.dataobject.TreeList.get_from_stream`, :meth:`~dendropy.dataobject.TreeList.get_from_path`, or :meth:`~dendropy.dataobject.TreeList.get_from_string`) to simulataneously instantiate and populate a :class:`~dendropy.dataobject.TreeList`.