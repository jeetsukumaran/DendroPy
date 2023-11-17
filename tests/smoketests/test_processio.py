from dendropy.utility.processio import Session
import io
import tempfile

from . import marksmoke as pytestmark

def test_start():
    read = Session()

    res = read.start("true")
    assert res is None
