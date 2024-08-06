import dendropy
from dendropy.calculate.profiledistance import MeasurementProfile, TreeProfile

from . import marksmoke as pytestmark


def test_add():
    prof = MeasurementProfile()

    res = prof.add(1)
    assert res is None


def test_set_data():
    prof = MeasurementProfile()

    res = prof.set_data([])
    assert res is None


def test_euclidean_distance():
    prof = MeasurementProfile()

    res = prof._euclidean_distance([], [])
    assert isinstance(res, float)


def test_distance():
    prof1 = MeasurementProfile()
    prof2 = MeasurementProfile()
    prof1.add(1)
    prof2.add(1)
    prof2.add(1)

    res = prof1.distance(prof2, 2)
    assert isinstance(res, float)


def test_compile():
    tree = dendropy.Tree()
    tree.seed_node.new_child(edge_length=1)
    prof = TreeProfile(tree)

    res = prof.compile(tree)
    assert res is None
