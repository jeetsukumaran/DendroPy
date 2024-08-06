import tempfile
from types import GeneratorType

from dendropy.utility.container import (
    CaseInsensitiveDict,
    DataTable,
    NormalizedBitmaskDict,
    OrderedCaselessDict,
    OrderedSet,
)

from . import marksmoke as pytestmark


def test_remove():
    set = OrderedSet([1])

    res = set.remove(1)
    assert res is None


def test_index():
    set = OrderedSet([1])

    res = set.index(1)
    assert isinstance(res, int)


def test_clear():
    set = OrderedSet([1])

    res = set.clear()
    assert res is None

def test_pop():
    set = NormalizedBitmaskDict()
    res = set.pop(1)
    assert isinstance(res, int)

def test_get():
    set = NormalizedBitmaskDict(fill_bitmask=0b11101001)
    res = set.get(1)
    assert res is None


def test_lower_items():
    dict = CaseInsensitiveDict()

    res = dict.lower_items()
    assert isinstance(res, GeneratorType)


def test_copy():
    dict = CaseInsensitiveDict()

    res = dict.copy()
    assert isinstance(res, CaseInsensitiveDict)


def test_caseless_copy():
    dict = OrderedCaselessDict()

    res = dict.copy()
    assert isinstance(res, OrderedCaselessDict)


def test_itervalues():
    dict = OrderedCaselessDict()

    res = dict.itervalues()
    assert isinstance(res, GeneratorType)


def test_iteritems():
    dict = OrderedCaselessDict()

    res = dict.iteritems()
    assert isinstance(res, GeneratorType)


def test_values():
    dict = OrderedCaselessDict()

    res = dict.values()
    assert isinstance(res, list)


def test_pop():
    dict = OrderedCaselessDict()

    res = dict.pop("A")
    assert res is None


def test_pop_item():
    dict = OrderedCaselessDict({"A": "A"})

    res = dict.popitem()
    assert isinstance(res, tuple)


def test_caseless_keys():
    dict = OrderedCaselessDict()

    res = dict.caseless_keys()
    assert isinstance(res, list)


def test_index():
    dict = OrderedCaselessDict({"A": "A"})

    res = dict.index("A")
    assert isinstance(res, int)


def test_keys():
    dict = OrderedCaselessDict()

    res = dict.keys()
    assert isinstance(res, list)


def test_caseless_clear():
    dict = OrderedCaselessDict()

    res = dict.clear()
    assert res is None


def test_has_key():
    dict = OrderedCaselessDict()

    res = dict.has_key("A")
    assert isinstance(res, bool)


def test_caseless_get():
    dict = OrderedCaselessDict()

    res = dict.get("A")
    assert res is None


def test_setdefault():
    dict = OrderedCaselessDict()

    res = dict.setdefault("A")
    assert res is None


def test_update():
    dict1 = OrderedCaselessDict()
    dict2 = OrderedCaselessDict()

    res = dict1.update(dict2)
    assert res is None


def test_fromkeys():
    dict = OrderedCaselessDict()

    res = dict.fromkeys([])
    assert isinstance(res, OrderedCaselessDict)


def test_row_value_iter():
    table = DataTable()

    res = table.row_value_iter("A")
    assert isinstance(res, GeneratorType)


def test_column_value_iter():
    table = DataTable()

    res = table.column_value_iter("A")
    assert isinstance(res, GeneratorType)


def test_write_csv():
    table = DataTable()

    with tempfile.NamedTemporaryFile() as file:
        res = table.write_csv(file.name)
    assert res is None


def test_num_columns():
    table = DataTable()

    res = table.num_columns()
    assert isinstance(res, int)
