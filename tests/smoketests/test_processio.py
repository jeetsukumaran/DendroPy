import tempfile

from dendropy.utility.processio import Session, SessionReader

from . import marksmoke as pytestmark

# TODO: https://github.com/jeetsukumaran/DendroPy/issues/179
# def test_enqueue_stream():
#     read = SessionReader("")
#     read.enqueue_stream()


def test_read():
    read = SessionReader("")

    res = read.read()
    assert res is None


def test_start():
    read = Session()

    res = read.start("true")
    assert res is None
