import tempfile

from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
from dendropy.utility.container import DataTable

from . import marksmoke as pytestmark


def test_as_data_table():
    mat = PhylogeneticDistanceMatrix()

    res = mat.as_data_table()
    assert isinstance(res, DataTable)


def test_write_csv():
    mat = PhylogeneticDistanceMatrix()

    with tempfile.NamedTemporaryFile() as file:
        res = mat.write_csv(file.name)
        assert res is None


def test_distances():
    mat = PhylogeneticDistanceMatrix()

    res = mat.distances()
    assert isinstance(res, list)
