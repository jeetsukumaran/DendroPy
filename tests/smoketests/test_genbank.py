from dendropy.interop.genbank import GenBankResourceStore

from . import pytestmark


def test_parsexml():
    s = GenBankResourceStore(1)

    res = s.parse_xml(string="<a></a>")
    assert isinstance(res, list)


def test_prepareids():
    s = GenBankResourceStore(1)

    res = s.prepare_ids([1])
    assert isinstance(res, list)
