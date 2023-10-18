import dendropy

# force tree to be read as rooted
mle_rooted = dendropy.Tree.get(
        path='pythonidae.mle.nex',
        schema='nexus',
        rooting='force-rooted')

# force tree to be read as unrooted
mle_rooted = dendropy.Tree.get(
        path='pythonidae.mle.nex',
        schema='nexus',
        rooting='force-unrooted')
