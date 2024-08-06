from types import GeneratorType

from dendropy.dataio.nexusprocessing import get_rooting_argument, group_ranges

from . import marksmoke as pytestmark


def test_group_ranges():
    res = group_ranges([[]])
    assert isinstance(res, GeneratorType)


def test_get_rooting_argument():
    res = get_rooting_argument(is_rooted=True)
    assert isinstance(res, str)
