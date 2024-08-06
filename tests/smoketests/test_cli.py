import tempfile

from dendropy.utility.cli import (
    compose_citation_for_program,
    confirm_overwrite,
    show_splash,
)

from . import marksmoke as pytestmark


def test_confirm_overwrite():
    with tempfile.NamedTemporaryFile() as file:
        res = confirm_overwrite(file.name, replace_without_asking=True)
    assert res is True


def test_show_splash():
    res = show_splash("a", "a", "a", "a", "a")
    assert res is None


def test_compose_citation_for_program():
    res = compose_citation_for_program("a", "a", "a", "a")
    assert isinstance(res, list)
