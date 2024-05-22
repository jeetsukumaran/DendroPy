from dendropy.utility.debug import dump_stack, get_calling_code_info

from . import marksmoke as pytestmark


def test_get_calling_code_info():
    res = get_calling_code_info(1)
    assert isinstance(res, tuple)


def test_dump_stack():
    res = dump_stack()
    assert res is None
