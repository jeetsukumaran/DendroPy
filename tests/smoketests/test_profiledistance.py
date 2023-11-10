import dendropy
from dendropy.calculate.profiledistance import MeasurementProfile, TreeProfile

from . import pytestmark


def test_add():
    p = MeasurementProfile()

    res = p.add(1)
    assert res is None


def test_set_data():
    p = MeasurementProfile()

    res = p.set_data([])
    assert res is None


def test_euclidean_distance():
    p = MeasurementProfile()

    res = p._euclidean_distance([], [])
    assert isinstance(res, float)


def test_distance():
    p1 = MeasurementProfile()
    p2 = MeasurementProfile()
    p1.add(1)
    p2.add(1)
    p2.add(1)

    res = p1.distance(p2, 2)
    assert isinstance(res, float)


def test_compile():
    t = dendropy.Tree()
    t.seed_node.new_child(edge_length=1)
    p = TreeProfile(t)

    res = p.compile(t)
    assert res is None
