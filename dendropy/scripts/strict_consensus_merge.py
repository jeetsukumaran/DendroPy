#!/usr/bin/env python
import sys
import copy
from dendropy import dataio
from dendropy import get_logger

_LOG = get_logger('scripts.strict_consensus_merge')
verbose = False
showTrees = False
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
        return ' -> '.join([str(splitToList(i.split, -1, True)) for i in nd])
    def collapse(self, leafIntersection):
        upperNd = self.upperNd
        if verbose:
            _LOG.debug('in collapse upperNd =')
            #upperNd.writeNewick(sys.stdout)
            _LOG.debug('\n') 
        higherPaths, collapseIrr = partitionChildrenWithLeafSet(upperNd, leafIntersection)
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
                    _LOG.debug('  %s' % ' '.join(map(str,splitToList(c.split, -1, True))))
            if upperNd.isInternal():
                upperNd.collapseEdge()
            if verbose:
                _LOG.debug( 'after collapse')
                for c in unp.children:
                    _LOG.debug('  %s' % ' '.join(map(str,splitToList(c.split, -1, True))))
        self.lowerNd._addChildren(collapseIrr)
        return higherPaths, collapseIrr
        
    def scmMergePath(self, toModifyNode, toDesPaths, toRecurse, leafIntersection):
        found = False
        collapseTM = False
        while not (found or collapseTM) and len(toDesPaths):
            collapseTD = False
            toIterate = copy.copy(toDesPaths)
            for desPathNumber, toDesPath in enumerate(toIterate):
                if verbose:
                    _LOG.debug('toDesPath#%d = %s' % (desPathNumber, str(toDesPath)))
                res = classifyPaths(self, toDesPath, leafIntersection)
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
                    higherPaths, collapsed  = toDesPath.collapse(leafIntersection)
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
            higherPaths, collapsed = self.collapse(leafIntersection) 
            return higherPaths
        return []
            
        

def partitionChildrenWithLeafSet(nd, leafIntersection):
    '''returns VerticalPath list for children that have leaves in the leafIntersection, and list of other children as a tuple.  
    Assumes 'nd' is a significant node in terms of the commonLeafSet'''
    lowerPattern = nd.split & leafIntersection
    assert(lowerPattern)
    paths = []
    irrelevantChildren = []
    for c in iterChildren(nd):
        if c.split & leafIntersection:
            assert(c.split & leafIntersection != lowerPattern)
            paths.append(VerticalPath(leafIntersection, nd, c))
        else:
            irrelevantChildren.append(c)
    return paths, irrelevantChildren
    
def moveInterveningNodesToRefEdge(tree, leafIntersection):
    '''Moves "irrelevant" nodes to the left side of the root maintaining the unrooted topology and root node object.
    Used to guarantee that the root is "significant" wrt leafIntersection and the path to the reference Node contains all intervening nodes'''
    root = tree.root
    keepMovingSplitRep = leafIntersection ^ root.lChild.split # if a node has this split, then there aren't any other children that have any of the intersection taxa
    while True:
        for c in iterChildren(root):
            if c.split & leafIntersection == keepMovingSplitRep: # slide root so that root.children are root.lChild.children
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
    
def rootPortionSCM(toModify, toDestroy, leafIntersection):
    '''Roots trees at a common leaf, and deals with collisions along that terminal path will be the root's lChild after calling the function.'''
    lowestCommonLeaf = getFirstBitAsIndex(leafIntersection)
    toModRefNode = toModify.attachAtLeafByIndex(lowestCommonLeaf)
    toDestroyRefNode = toDestroy.attachAtLeafByIndex(lowestCommonLeaf)
    moveInterveningNodesToRefEdge(toModify, leafIntersection)
    moveInterveningNodesToRefEdge(toDestroy, leafIntersection)
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
    
def classifyPaths(firPath, secPath, leafIntersection):
    #print 'classifyPaths(' +  str(splitToList(firPath.split, -1, True)) + ', ' +  str(splitToList(secPath.split, -1, True)) + ')'
    comp = compareSplits(firPath.split, secPath.split, leafIntersection, False)[0]
    if comp != SplitCompatEnum.compat:
        return comp
    intersec = firPath.split & secPath.split
    if intersec == firPath.split:
        return SplitCompatEnum._firstIsSubset
    return intersec == secPath.split and SplitCompatEnum._secondIsSubset or SplitCompatEnum.compat
    
def scmMergeSubTree(toModifyNode, toDestroyNode, leafIntersection):
    #print 'leafIntersection =', splitToList(leafIntersection, -1, True)
    toModPaths, toModOtherNds = partitionChildrenWithLeafSet(toModifyNode, leafIntersection)
    toDesPaths, toDestroyOtherNds = partitionChildrenWithLeafSet(toDestroyNode, leafIntersection)
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
            newPaths.extend(toModPath.scmMergePath(toModifyNode, toDesPaths, toRecurse, leafIntersection))
            if verbose:
                _LOG.debug('\n')
        if not newPaths:
            break
        else:
            toModPaths = newPaths
    for p in toRecurse:
        scmMergeSubTree(p[0], p[1], leafIntersection)
    
    
def add_to_scm(treeToModify, treeToDestroy):
    leafIntersection = treeToModify.root.split & treeToDestroy.root.split
    _LOG.debug('leaf intersection = %s' % ' '.join(map(str,splitToList(leafIntersection, -1, True))))
    nCommonLeaves = countBits(leafIntersection)
    if nCommonLeaves < 2:
        _LOG.error('trees must have at least 2 common leaves')
        raise ValueError, 'trees must have at least 2 common leaves'
    if nCommonLeaves == 2: # return a polytomy
        treeToModify.root.collapseClade()
        treeToDestroy.root.collapseClade()
        treeToModify.root._stealChildren([c for c in iterChildren(treeToDestroy.root) if not (leafIntersection & c.split)])
    else:
        rootPortionSCM(treeToModify, treeToDestroy, leafIntersection)
        if verbose and showTrees:
            treeToModify.show('toModPostRoot.dot')
            treeToDestroy.show('toDestroyPostRoot.dot')
        if verbose:
            _LOG.debug('treeToDestroy.root.children = ')
            for c in treeToDestroy.root.children:
                _LOG.debug('  %s'  %  ' '.join(map(str,splitToList(c.split, -1, True))))
            
        # prune the reference edge that we have already dealt with
        toReattach = treeToModify.root.lChild
        toReattach.pruneSelf(collapseDegTwo = False)
        treeToDestroy.root.lChild.pruneSelf(collapseDegTwo = False)
        scmMergeSubTree(treeToModify.root, treeToDestroy.root, leafIntersection)
            # reattach the root's left child
        toReattach.rSib = treeToModify.root.lChild
        toReattach.par = treeToModify.root
        treeToModify.root.lChild = toReattach
        del toReattach
    treeToModify.root.refreshSplits()
    treeToModify.taxaManager |= treeToDestroy.taxaManager
    if hasattr(treeToModify, 'leafSet'):
        del treeToModify.leafSet
    return treeToModify



def strict_consensus_merger(self, tree_list, copy_trees=True):
    """Returns a tree that is the strict consensus merger of the input trees.
    
    If copy_trees is True then the trees will be copied before the merger 
    operation, if the `copy_trees` is False then the input trees will be 
    destroyed by the operation (and the modified first tree will be returned).
    """
    if copy_tree:
        tree_list = [copy.copy(i) for i in tree_list]
    nTrees = len(tree_list)
    _LOG.debug('%d Trees to merge:\n%s\n' % (nTrees, '\n'.join([i.m_newick for i in tree_list])))
    if nTrees < 2:
        return tree_list[0]
    tree_iter = iter(tree_list)
    target_tree = tree_iter.next()
    encode_splits(target_tree)
    for next_tree in tree_iter:
        encode_splits(next_tree)
        add_to_scm(target_tree, next_tree)
    return target_tree


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


