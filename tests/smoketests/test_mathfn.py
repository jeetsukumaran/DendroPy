from dendropy.calculate.mathfn import LCM, lcm

from . import marksmoke as pytestmark


def test_lcm():
    res = lcm(32, 67)
    assert isinstance(res, int)


def test_LCM():
    res = LCM(1, 2, 3, 4, 5)
    assert isinstance(res, int)
