from dendropy.interop.genbank import GenBankResourceStore

from . import marksmoke as pytestmark


def test_parsexml():
    store = GenBankResourceStore(1)

    res = store.parse_xml(string="<a></a>")
    assert isinstance(res, list)


def test_prepareids():
    store = GenBankResourceStore(1)

    res = store.prepare_ids([1])
    assert isinstance(res, list)
