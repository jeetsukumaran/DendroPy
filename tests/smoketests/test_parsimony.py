import dendropy
from dendropy.model.parsimony import fitch_up_pass

from . import pytestmark


def test_fitch_up_pass():
    t = dendropy.Tree()

    res = fitch_up_pass(t.preorder_node_iter())
    assert res is None
