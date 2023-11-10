from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
from dendropy.utility.container import DataTable

from . import pytestmark


def test_as_data_table():
    m = PhylogeneticDistanceMatrix()

    res = m.as_data_table()
    assert isinstance(res, DataTable)


def test_write_csv():
    m = PhylogeneticDistanceMatrix()

    # res = m.write_csv("a.out")
    # assert res == None


def test_distances():
    m = PhylogeneticDistanceMatrix()

    res = m.distances()
    assert isinstance(res, list)
