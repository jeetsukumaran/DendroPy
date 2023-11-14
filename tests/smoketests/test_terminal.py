from dendropy.utility.terminal import ttysize, win_terminal_width

from . import marksmoke as pytestmark


def test_ttysize():
    res = ttysize()
    assert isinstance(res, tuple)
