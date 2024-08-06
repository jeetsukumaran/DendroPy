import io
import tempfile

from dendropy.utility.filesys import (
    LineReadingThread,
    find_executable,
    find_files,
    glob_match,
)

from . import marksmoke as pytestmark


def test_wait_for_file_to_appear():
    thread = LineReadingThread()

    with tempfile.NamedTemporaryFile() as file:
        res = thread.wait_for_file_to_appear(file.name)
    assert isinstance(res, bool)


def test_open_file_when_exists():
    thread = LineReadingThread()

    with tempfile.NamedTemporaryFile() as file:
        res = thread.open_file_when_exists(file.name)
    assert isinstance(res, io.IOBase)


def test_run():
    thread = LineReadingThread()
    res = thread.run()
    assert res is None


def test_keep_going():
    thread = LineReadingThread()
    res = thread.keep_going("A")
    assert isinstance(res, bool)


def test_glob_match():
    res = glob_match("A", "A")
    assert isinstance(res, bool)


def test_find_files():
    res = find_files("A")
    assert isinstance(res, list)


def test_find_executable():
    res = find_executable("A")
    assert res is None
