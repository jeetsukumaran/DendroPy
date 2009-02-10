#!/usr/bin/env python
import sys
import copy
from dendropy import dataio
from dendropy.splits import encode_splits, split_to_list, count_bits, lowest_bit_only
from dendropy import get_logger
from dendropy.treemanip import collapse_clade

_LOG = get_logger('scripts.strict_consensus_merge')
verbose = False


from dendropy.utils import NormalizedBitmaskDict

def reroot_on_lowest_common_index_path(t, common_mask):
    """This operation is only for unrooted trees that are being merged using
    SCM. The path the separates the lowest index taxon in the leaf set 
    intersection is placed as the first descendant path of the "seed_node" for
    the tree.
    This assures that all split representations are oriented in the same way 
    for subsequent operations.
    The mask most contain more that 2 bits (there must be an internal node in
    the tree that is has degree > 2  wrt the common leafset).
    """
    l = lowest_bit_only(common_mask)
    assert(l > 0)
    assert(count_bits(common_mask) > 2)
    # start at the lowest leaf in common.
    curr_n = t.split_edges[l].head_node
    # walk back toward the root until we find a node that has another bit
    p = curr_n.parent_node
    while p:
        if (p.edge.clade_mask & common_mask) != l:
            break
        curr_n = p
        p = curr_n.parent_node
    if (curr_n.edge.clade_mask & common_mask) == l:
        # we did not make it to the root.  Make curr_n, the first_child of the root
        t.to_outgroup_position(curr_n, flip_splits=True, suppress_deg_two=True)
        return
    # we hit the root, now we walk up the tree, to find the a relevant internal
    lowest_on_path_to_l = curr_n
    taxa_mask = t.seed_node.clade_mask
    comp_mask = (~common_mask) & taxa_mask
    children = curr_n.child_nodes()
    assert(len(children) > 1)
    nd_source = iter(children)
    try:
        while True:
            c = nd_source.next()
            cm = c.edge.clade_mask
            masked_cm = cm & comp_mask
            if masked_cm:
                if masked_cm == comp_mask:
                    curr_n = c
                    children = curr_n.child_nodes()
                    assert(len(children) > 1)
                    nd_source = iter(children)
                else:
                    break
    except StopIteration:
        raise AssertionError("Reaching supposedly unreachable code")
    
    if curr_n is not t.seed_node:
        # We have found the first relevant internal node, we want to make it 
        #   the root.  We can do this by putting one of its children into the
        #   "outgroup position" and then putting the path to lowest commond
        #   leaf in the outgroup position (this last operation is just a 
        #   rearrangement of the order of children in the root.
        children = curr_n.child_nodes()
        assert(len(children) > 1)
        p = curr_n.parent
        t.to_outgroup_position(children[0], flip_splits=True, suppress_deg_two=True)
        t.to_outgroup_position(p, flip_splits=True, suppress_deg_two=True)
    else:
        # if the root first relevant, node then we just make the path leading
        #   to the lowest index node the first child of the root
        t.to_outgroup_position(lowest_on_path_to_l, flip_splits=True, suppress_deg_two=True)



def add_to_scm(to_modify, to_consume, rooted=False):
    """Adds the tree `to_consume` to the tree `to_modify` in a strict consensus
    merge operation.  Both trees must have had encode_splits called on them."""
    if rooted:
        raise NotImplementedError("rooted form of add_to_scm not implemented")
    to_mod_root = to_modify.seed_node
    to_mod_split = to_mod_root.edge.clade_mask

    to_consume_root = to_consume.seed_node
    to_consume_split = to_consume_root.edge.clade_mask

    leaf_intersection = to_mod_split & to_consume_split
    n_common_leaves = count_bits(leaf_intersection)
    if n_common_leaves < 2:
        _LOG.error('trees must have at least 2 common leaves')
        raise ValueError('trees must have at least 2 common leaves')
    if n_common_leaves == 2: 
        # SCM with 2 leaves in common results in a polytomy
        collapse_clade(to_mod_root)
        collapse_clade(to_consume_root)
        leaves_to_steal = [c for c in to_consume_root.child_nodes() if not (leaf_intersection & c.edge.clade_mask)]
        for leaf in leaves_to_steal:
            to_mod_root.add_child(leaf)
            to_mod_root.edge.clade_mask |= leaf.edge.clade_mask
        return
    to_mod_relevant_splits = {}
    to_consume_relevant_splits = {}
    if not rooted:
        reroot_on_lowest_common_index_path(to_modify, leaf_intersection)
        reroot_on_lowest_common_index_path(to_consume, leaf_intersection)
        to_mod_root = to_modify.seed_node
        assert(to_mod_root.edge.clade_mask == to_mod_split)
        to_consume_root = to_consume.seed_node
        assert(to_consume_root.edge.clade_mask == to_consume_split)

    for s, e in to_modify.split_edges:
        masked = s & leaf_intersection
        e_list = to_mod_relevant_splits.setdefault(masked, [])
        e_list.append(e)
    for s, e in to_consume_relevant_splits.split_edges:
        masked = s & leaf_intersection
        e_list = to_consume_relevant_splits.setdefault(masked, [])
        e_list.append(e)
    
    to_steal = [i for i in to_consume_root.children() if (i.edge.clade_mask & leaf_intersection) == 0]
    for child in to_steal:
        to_mod_root.add_child(child)
        to_mod_root.edge.clade_mask |= child.edge.clade_mask
    for split, path in to_consume_relevant_splits.iteritems():
        code=here
        
        
    
def strict_consensus_merge(tree_list, taxa_block, copy_trees=False, rooted=False):
    """Returns a tree that is the strict consensus merger of the input trees.
    
    If copy_trees is True then the trees will be copied before the merger 
    operation, if the `copy_trees` is False then the input trees will be 
    destroyed by the operation (and the modified first tree will be returned).
    """
    if copy_trees:
        tree_list = [copy.copy(i) for i in tree_list]
    nTrees = len(tree_list)
    _LOG.debug('%d Trees to merge:\n%s\n' % (nTrees, '\n'.join([str(i) for i in tree_list])))
    if nTrees < 2:
        return tree_list[0]
    tree_iter = iter(tree_list)
    to_modify = tree_iter.next()
    encode_splits(to_modify, taxa_block=taxa_block)

    if rooted:
        raise NotImplementedError("Rooted SCM is not implemented")
    for to_consume in tree_iter:
        encode_splits(to_consume, taxa_block=taxa_block)
        add_to_scm(to_modify, to_consume)

    return to_modify


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-p", "--prec", dest="prec", default=0.00001, 
        type="float",
        help="The precision of the comparison that the node depth is equal regardless of tip to node path.")
    (options, args) = parser.parse_args()
    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified")
    if len(args) == 1:
        newick = args[0]
    else:
        newick = sys.stdin.readline()

    prec = options.prec
    tree = dataio.trees_from_string(string=newick, format="NEWICK")[0]
    dendropy.trees.add_depth_to_nodes(tree, options.prec)

    sys.stdout.write("%f\n" % tree.seed_node.depth)




x = """showTrees = False
conflictDestroysAllSubTreeStructure = False

class VerticalPath:
    def _findUpperNode(nd, commonLeafSet):
        lowerPattern = nd.split & commonLeafSet
        intervening = []
        while True:
            if nd.isLeaf():
                return nd, intervening
            nextNd = None
            for c in iterChildren(nd):
                inters = c.split & commonLeafSet
                if inters:
                    if inters == lowerPattern:
                        nextNd = c
                    else:
                        return nd, intervening
            intervening.append(nd)
            assert(nextNd)
            nd = nextNd
    _findUpperNode = staticmethod(_findUpperNode)
    
    '''A "vertical" path through a tree containing an arbitrary number of intervening nodes (in order of root end to tip).'''
    def __init__(self, mask , lowerNd, lowerNdChild):
        self.lowerNd = lowerNd
        upperNd, trunkNdList = VerticalPath._findUpperNode(lowerNdChild, mask)
        self.upperNd = upperNd
        self.trunkNodes = trunkNdList
        if mask and upperNd:
            self.split = mask & (upperNd.split)
    def __str__(self):
        nd = [self.lowerNd] + self.trunkNodes + [self.upperNd]
        return ' -> '.join([str(split_to_list(i.split, -1, True)) for i in nd])

    def collapse(self, leaf_intersection):
        upperNd = self.upperNd
        if verbose:
            _LOG.debug('in collapse upperNd =')
            #upperNd.writeNewick(sys.stdout)
            _LOG.debug('\n') 
        higherPaths, collapseIrr = partitionChildrenWithLeafSet(upperNd, leaf_intersection)
        #if verbose:
            #_LOG.debug('higherPaths = ')
            #for h in higherPaths:
            #    h.upperNd.writeNewick(sys.stdout)
            #     
        for h in higherPaths:
            h.lowerNd = self.lowerNd
        for c in collapseIrr:
            c.pruneSelf()
        if len(self.trunkNodes) > 0:
            upperNd.pruneSelf(False)
            lowestTrunk = self.trunkNodes[0]
            lowestTrunk.pruneSelf(False)
            if conflictDestroysAllSubTreeStructure:
                lowestTrunk.collapseClade()
            else:
                for n in self.trunkNodes[1:]:
                    n.collapseEdge()
            collapseIrr.extend(lowestTrunk.children)
            for h in higherPaths:
                self.lowerNd._addChild(len(h.trunkNodes) and h.trunkNodes[0] or h.upperNd)
        else:
            if verbose:
                _LOG.debug('upperNd.par.children = ')
                unp = upperNd.par
                for c in unp.children:
                    _LOG.debug('  %s' % ' '.join(map(str,split_to_list(c.split, -1, True))))
            if upperNd.isInternal():
                upperNd.collapseEdge()
            if verbose:
                _LOG.debug( 'after collapse')
                for c in unp.children:
                    _LOG.debug('  %s' % ' '.join(map(str,split_to_list(c.split, -1, True))))
        self.lowerNd._addChildren(collapseIrr)
        return higherPaths, collapseIrr
        
    def scmMergePath(self, toModifyNode, toDesPaths, toRecurse, leaf_intersection):
        found = False
        collapseTM = False
        while not (found or collapseTM) and len(toDesPaths):
            collapseTD = False
            toIterate = copy.copy(toDesPaths)
            for desPathNumber, toDesPath in enumerate(toIterate):
                if verbose:
                    _LOG.debug('toDesPath#%d = %s' % (desPathNumber, str(toDesPath)))
                res = classifyPaths(self, toDesPath, leaf_intersection)
                if res == SplitCompatEnum.identical:
                    _LOG.debug('identical')
                    if len(toDesPath.trunkNodes) > 0:
                        if len(self.trunkNodes) > 0: # collision
                            _LOG.debug('collision')
                            self.upperNd.pruneSelf(False)
                            if conflictDestroysAllSubTreeStructure:
                                self.trunkNodes[0].collapseDescendants()
                            else:
                                for n in self.trunkNodes[1:]:
                                    n.collapseEdge()
                            self.trunkNodes[0]._addChild(self.upperNd)
                            toDesPath.upperNd.pruneSelf(collapseDegTwo = False)
                            if conflictDestroysAllSubTreeStructure:
                                toDesPath.trunkNodes[0].collapseDescendants()
                            else:
                                for n in toDesPath.trunkNodes[1:]:
                                    n.collapseEdge()
                            self.trunkNodes[0]._addChildren(toDesPath.trunkNodes[0].children)
                            toDesPath.lChild = None
                            toDesPath.trunkNodes[0].pruneSelf()
                        else:
                            self.upperNd._swapPlaces(toDesPath.trunkNodes[0])
                            toDesPath.upperNd._swapPlaces(self.upperNd)
                    if self.upperNd.isInternal():
                        toRecurse.append((self.upperNd, toDesPath.upperNd))
                    toDesPaths.remove(toDesPath)
                    found = True
                    break
                elif res == SplitCompatEnum.incompat:
                    _LOG.debug('incompat')
                    collapseTM = True
                    collapseTD = True
                elif res == SplitCompatEnum._firstIsSubset:
                    _LOG.debug('_firstIsSubset')
                    collapseTD = True
                elif res == SplitCompatEnum._secondIsSubset:
                    _LOG.debug('_secondIsSubset')
                    collapseTM = True
                if collapseTD:
                    higherPaths, collapsed  = toDesPath.collapse(leaf_intersection)
                    toModifyNode._stealChildren(collapsed)
                    toDesPaths.remove(toDesPath)
                    toDesPaths.extend(higherPaths)
                    break
                if collapseTM:
                    break
        if not found:
            _LOG.debug('not found')
            collapseTM = True
        if collapseTM:
            higherPaths, collapsed = self.collapse(leaf_intersection) 
            return higherPaths
        return []
            
        

def partitionChildrenWithLeafSet(nd, leaf_intersection):
    '''returns VerticalPath list for children that have leaves in the leaf_intersection, and list of other children as a tuple.  
    Assumes 'nd' is a significant node in terms of the commonLeafSet'''
    lowerPattern = nd.split & leaf_intersection
    assert(lowerPattern)
    paths = []
    irrelevantChildren = []
    for c in iterChildren(nd):
        if c.split & leaf_intersection:
            assert(c.split & leaf_intersection != lowerPattern)
            paths.append(VerticalPath(leaf_intersection, nd, c))
        else:
            irrelevantChildren.append(c)
    return paths, irrelevantChildren
    
def moveInterveningNodesToRefEdge(tree, leaf_intersection):
    '''Moves "irrelevant" nodes to the left side of the root maintaining the unrooted topology and root node object.
    Used to guarantee that the root is "significant" wrt leaf_intersection and the path to the reference Node contains all intervening nodes'''
    root = tree.root
    keepMovingSplitRep = leaf_intersection ^ root.lChild.split # if a node has this split, then there aren't any other children that have any of the intersection taxa
    while True:
        for c in iterChildren(root):
            if c.split & leaf_intersection == keepMovingSplitRep: # slide root so that root.children are root.lChild.children
                c.pruneSelf(collapseDegTwo = False)
                nextRootChildren = c.children
                c.lChild = None
                c._addChildren(root.children)
                root.lChild = None
                c.split = root.split ^ c.split
                root._addChild(c)
                root._addChildren(nextRootChildren)
                break
            if c.split & keepMovingSplitRep:
                return
    
def _root_portion_SCM(result_tree, tree_to_merge, leaf_intersection):
    '''roots trees at a common leaf, and deals with collisions along that 
    terminal path. 
    
    The leaf with the lowest index will be the root's first child after calling
    the function.'''
    lowest_leaf_mask = lowest_bit_only(leaf_intersection)
    lowest_leaf_in_result_tree = result_tree.split_edges[lowest_leaf_mask].head_node
    result_tree.reroot_at(lowest_leaf_in_result_tree.parent_node, flip_splits=True)

    lowest_leaf_in_to_merge = tree_to_merge.split_edges[lowest_leaf_mask].head_node
    tree_to_merge.reroot_at(lowest_leaf_in_to_merge.parent_node, flip_splits=True)

    toModRefNode = toModify.attachAtLeafByIndex(lowestCommonLeaf)
    toDestroyRefNode = toDestroy.attachAtLeafByIndex(lowestCommonLeaf)
    moveInterveningNodesToRefEdge(toModify, leaf_intersection)
    moveInterveningNodesToRefEdge(toDestroy, leaf_intersection)
    if toModRefNode.par.par:
        if toDestroyRefNode.par.par: #collision
            if conflictDestroysAllSubTreeStructure:
                toDestroy.root.lChild.collapseDescendants()
                toModify.root.lChild.collapseDescendants()
            else: #reduce to two edges on the refNode to root path
                while toDestroyRefNode.par.par.par:
                    toDestroyRefNode.par.collapseEdge()
                while toModRefNode.par.par.par:
                    toModRefNode.par.collapseEdge()
            p = toModRefNode.par
            toSteal = [c for c in iterChildren(toDestroyRefNode.par) if c != toDestroyRefNode]
            toDestroyRefNode.par.lChild = toDestroyRefNode
            toDestroyRefNode.rSib = None
            p._addChildren(toSteal)
    else:
        if toDestroyRefNode.par.par:
            toModRefNode._swapPlaces(toDestroy.root.lChild)
            toDestroyRefNode._swapPlaces(toModRefNode) # for completeness, we'll move the "correct" reference node back onto it's original tree
    
def classifyPaths(firPath, secPath, leaf_intersection):
    #print 'classifyPaths(' +  str(split_to_list(firPath.split, -1, True)) + ', ' +  str(split_to_list(secPath.split, -1, True)) + ')'
    comp = compareSplits(firPath.split, secPath.split, leaf_intersection, False)[0]
    if comp != SplitCompatEnum.compat:
        return comp
    intersec = firPath.split & secPath.split
    if intersec == firPath.split:
        return SplitCompatEnum._firstIsSubset
    return intersec == secPath.split and SplitCompatEnum._secondIsSubset or SplitCompatEnum.compat
    
def scmMergeSubTree(toModifyNode, toDestroyNode, leaf_intersection):
    #print 'leaf_intersection =', split_to_list(leaf_intersection, -1, True)
    toModPaths, toModOtherNds = partitionChildrenWithLeafSet(toModifyNode, leaf_intersection)
    toDesPaths, toDestroyOtherNds = partitionChildrenWithLeafSet(toDestroyNode, leaf_intersection)
    if len(toDestroyOtherNds) > 0:
        toModifyNode._stealChildren(toDestroyOtherNds)
    toRecurse = []
    #print 'toModPaths =\n\t', '\n\t'.join([str(p) for p in toModPaths])
    while 1:
        newPaths = []
        for toModPath in toModPaths:
            if verbose:
                _LOG.debug('toModPath = %s' % toModPath)
                _LOG.debug('toDesPaths = %s' % '\n     '.join(map(str, toDesPaths)))
            newPaths.extend(toModPath.scmMergePath(toModifyNode, toDesPaths, toRecurse, leaf_intersection))
            if verbose:
                _LOG.debug('\n')
        if not newPaths:
            break
        else:
            toModPaths = newPaths
    for p in toRecurse:
        scmMergeSubTree(p[0], p[1], leaf_intersection)
    
    
def add_to_scm(tree_to_modify, tree_to_destroy, taxa_block):
    print str(tree_to_modify)
    print str(tree_to_destroy)
    print('leaf intersection = %s' % ' '.join([str(i) for i in split_to_list(leaf_intersection, -1, True)]))
    if n_common_leaves == 2: # return a polytomy
        collapse_clade(to_mod_root)
        collapse_clade(to_destroy_root)
        leaves_to_steal = [c for c in to_destroy_root.child_nodes() if not (leaf_intersection & c.edge.clade_mask)]
        for leaf in leaves_to_steal:
            to_mod_root.add_child(leaf)
    else:
        _root_portion_SCM(tree_to_modify, tree_to_destroy, leaf_intersection)
        if verbose and showTrees:
            tree_to_modify.show('toModPostRoot.dot')
            tree_to_destroy.show('toDestroyPostRoot.dot')
        if verbose:
            _LOG.debug('tree_to_destroy.root.children = ')
            for c in tree_to_destroy.root.children:
                _LOG.debug('  %s'  %  ' '.join(map(str,split_to_list(c.split, -1, True))))
            
        # prune the reference edge that we have already dealt with
        toReattach = tree_to_modify.root.lChild
        toReattach.pruneSelf(collapseDegTwo = False)
        tree_to_destroy.root.lChild.pruneSelf(collapseDegTwo = False)
        scmMergeSubTree(tree_to_modify.root, tree_to_destroy.root, leaf_intersection)
            # reattach the root's left child
        toReattach.rSib = tree_to_modify.root.lChild
        toReattach.par = tree_to_modify.root
        tree_to_modify.root.lChild = toReattach
        del toReattach
    encode_splits(tree_to_modify, taxa_block=taxa_block)
"""
