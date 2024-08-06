import tempfile

from dendropy.utility.bibtex import BibTexDb, BibTexEntry

from . import marksmoke as pytestmark

# TODO https://github.com/jeetsukumaran/DendroPy/issues/179
# def test_parse_text():
#     entry = BibTexEntry(citation="@asdfdsfadsf")
#     entry.parse_text()

# def test_fields_as_dict():
#     entry = BibTexEntry()
#     entry.fields_as_dict()


def test_add_from_file():
    db = BibTexDb()
    with tempfile.NamedTemporaryFile() as file:
        res = db.add_from_file(file.name)
    assert res is None


def test_add_from_text():
    db = BibTexDb()
    with tempfile.NamedTemporaryFile() as file:
        res = db.add_from_text("A")
    assert res is None
