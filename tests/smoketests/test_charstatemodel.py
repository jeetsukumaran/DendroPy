from dendropy.datamodel.charstatemodel import StateIdentity

from . import marksmoke as pytestmark


def test_is_exact_correspondence():
    id1 = StateIdentity()
    id2 = StateIdentity()

    res = id1.is_exact_correspondence(id2)
    assert isinstance(res, bool)
