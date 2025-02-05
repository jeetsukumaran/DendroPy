from dendropy.utility.textprocessing import (
    camel_case,
    format_dict_table,
    format_dict_table_rows,
    snake_case,
)

from . import marksmoke as pytestmark


def test_camel_case():
    res = camel_case("a")
    assert isinstance(res, str)


def test_snake_case():
    res = snake_case("a")
    assert isinstance(res, str)


def test_format_dict_table():
    res = format_dict_table("a")
    assert isinstance(res, str)


def test_format_dict_table_rows():
    res = format_dict_table_rows("a")
    assert isinstance(res, str)
