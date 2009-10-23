****************
Traversing Trees
****************

Trees can be traversed in pre-order, post-order, or level-order, over nodes or edges.

The following example demonstrates tree traversal. It calculates the ages of nodes (i.e., the node depths) and assigns the value to an attribute, ``age``, on each node. We traverse the tree in postorder, visiting children first. This way, for every node that we visit we are guaranteed that the child nodes already have their ages calculated, and so to get the age of the current node we just need to add the age of one of its child nodes to the edge connecting the current node to the child node. For this to be fully valid, the tree needs to ultrametric, which would mean that it would not matter which child node we picked.


.. topic:: Decorating Nodes with Node Ages   
    :class: code-recipe
    
    ::    
    
        def add_ages_to_nodes(tree, 
                              check_prec=0.0000001):
            """
            Takes an ultrametric `tree` and adds a attribute `age` to
            each node, with the value equal to the sum of edge lengths
            from the node to the tips. If the lengths of different paths
            to the node differ by more than `ultrametricity_prec`, then
            a ValueError exception will be raised indicating deviation
            from ultrametricity. If `ultrametricity_prec` is negative or
            False, then this check will be skipped.
            """
            node = None    
            for node in tree.postorder_node_iter():
                ch = node.child_nodes()
                if len(ch) == 0:
                    node.age = 0.0
                else:
                    first_child = ch[0]
                    node.age = first_child.age + first_child.edge.length)
                    if not (check_prec < 0 \
                            or check_prec == False):
                        for nnd in ch[1:]:
                            ocnd = nnd.age + nnd.edge.length
                            if abs(node.age - ocnd) > check_prec:
                                raise ValueError("Tree is not ultrametric")
            if node is None:
                raise ValueError("Empty tree encountered") 

The above example is actually based on a built-in method of the ``Tree`` class, ``add_ages_to_nodes()``, so in actual practice, if you do want to annotate nodes with their ages, you will not need to write such a function yourself, as you would use ``Tree.add_ages_to_nodes()`` directly. The following example shows how you might use this method to annotate nodes with their ages, and then report the list of ages as well as the age of the root.   

.. topic:: Decorating Nodes with Node Ages Using ``Tree.add_ages_to_nodes()`` 
    :class: code-recipe
    
    ::
                
        #! /usr/bin/env python
        
        from dendropy import datasets
        d = datasets.Dataset()
        d.read(open("results.tre", "rU"), "newick")
        for idx, tree in enumerate(d.trees_blocks[0]):
            tree.add_ages_to_nodes()
            node_ages = [node.age for node in tree.postorder_node_iter()]
            print("\nTree %d" % idx)
            print("    Node Ages: %s" % str(node_ages))
            print("    Age of Root: %f" % tree.seed_node.age)

                       
The following shows how you might calculate the total length of trees by visiting every edge and summing their lengths:

.. topic:: Calculating Tree Length
    :class: code-recipe
    
    ::

        #! /usr/bin/env python
        
        from dendropy import datasets
        
        def tree_length(tree):
            """Returns sum of branch lengths on tree."""
            total_length = 0
            for e in tree.postorder_edge_iter():
                if e.length is not None:
                    total_length += e.length
            return total_length
        
        d = datasets.Dataset()
        d.read( open("primates.tre", "ru"), "NEXUS" )
        for tb in d.trees_blocks:
            for t in tb:
                print("Tree Block '%s', Tree '%s': Length = %f" 
                        % (tb.label, t.label, tree_length(t)))

Because the ``length`` attribute of the root edge (i.e., the ``edge`` attribute of ``Tree.seed_node``) of an unrooted tree will be ``None``, we explicitly verify that each ``Edge`` object's ``length`` attribute is not ``None`` before adding to the sum.

The ``tree_length()`` function above could also be implemented by visiting nodes instead of edges::

    def tree_length(tree):
        """Returns sum of branch lengths on tree."""
        total_length = 0
        for n in tree.postorder_node_iter():
            if n.edge.length is not None:
                total_length += n.edge.length
        return total_length
