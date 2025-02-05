import dendropy
from dendropy.dataio.ioservice import DataReader
from dendropy.dataio.phylipreader import PhylipReader

from . import marksmoke as pytestmark


def test_describe_model():
    reader = PhylipReader()

    res = reader.describe_mode()
    assert isinstance(res, str)
