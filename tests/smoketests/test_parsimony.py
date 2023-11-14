import dendropy
from dendropy.model.parsimony import fitch_up_pass

from . import marksmoke as pytestmark


def test_fitch_up_pass():
    tree = dendropy.Tree()

    res = fitch_up_pass(tree.preorder_node_iter())
    assert res is None
