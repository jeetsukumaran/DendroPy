*******************************
Creating and Reading Tree Lists
*******************************

Multiple :class:`~dendropy.dataobject.Tree` objects referencing the same taxa can be aggregated into a  :class:`~dendropy.dataobject.TreeList`.
You can build up a :class:`~dendropy.dataobject.TreeList` object by instantiating an empty one, and individually appending  :class:`~dendropy.dataobject.Tree` objects to it::


