from dendropy.utility.processio import Session, SessionReader
import io
import tempfile

from . import marksmoke as pytestmark

# TODO: https://github.com/jeetsukumaran/DendroPy/issues/179
# def test_enqueue_stream():

#     with open("file") as file:
#         #file.readline()
#         read = SessionReader(file)
#         read.enqueue_stream()


def test_read():
    with open("file") as file:
        read = SessionReader(file)
        res = read.read()
        
    assert res is None


def test_start():
    read = Session()

    res = read.start("true")
    assert res is None
