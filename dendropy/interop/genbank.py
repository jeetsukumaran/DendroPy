#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Wrappers for retrieving data from GenBank.
"""

import sys
import re
import urllib
import dendropy
from dendropy.interop import entrez
from xml.etree import ElementTree
from dendropy.utility import container
from dendropy.utility import error

GENBANK_ANNOTATION_PREFIX = "genbank"
GENBANK_ANNOTATION_NAMESPACE = "http://www.ncbi.nlm.nih.gov/dtd/INSD_INSDSeq.mod.dtd"

##############################################################################
## GenBank Resources

class GenBankResourceStore(object):
    """
    Base GenBank data acquisition and manager class.
    """

    def parse_xml(**kwargs):
        if "stream" in kwargs and "string" in kwargs:
            raise TypeError("Cannot specify both 'stream' and 'string'")
        elif "stream" in kwargs:
            tree = ElementTree.parse(kwargs["stream"])
            root = tree.getroot()
        elif "string" in kwargs:
            # from StringIO import StringIO
            # s = StringIO(kwargs["string"])
            # try:
            #     root = ElementTree.parse(s)
            # except:
            #     print kwargs["string"]
            #     raise
            s = kwargs["string"]
            try:
                root = ElementTree.fromstring(s)
            except:
                sys.stderr.write(s)
                raise
        else:
            raise TypeError("Must specify exactly one of 'stream' or 'string'")
        gb_recs = []
        for seq_set in root.iter("INSDSet"):
            for seq in seq_set.iter("INSDSeq"):
                gb_rec = GenBankAccessionRecord(xml=seq)
                gb_recs.append(gb_rec)
        return gb_recs
    parse_xml = staticmethod(parse_xml)

    def fetch_xml(db, ids, prefix=None, email=None, as_stream=False):
        stream = entrez.efetch(db=db,
                ids=ids,
                rettype='gbc',
                retmode='xml',
                email=email)
        if as_stream:
            return stream
        else:
            return stream.read()
    fetch_xml = staticmethod(fetch_xml)

    def prepare_ids(ids, prefix=None):
        if prefix is None:
            prefix = ""
        # for idx, i in enumerate(ids):
        #     ids[idx] = "%s%s" % (prefix, i)
        ids = ["%s%s" % (prefix, i) for i in ids]
        return ids
    prepare_ids = staticmethod(prepare_ids)

    class AccessionFetchError(Exception):
        def __init__(self, missing, response=None):
            if response is None:
                response = ""
            else:
                response = "\n\nServer response was:\n\n%s" % response
            missing_desc = ", ".join([str(s) for s in missing])
            Exception.__init__(self, "\n\nFailed to retrieve accessions: %s%s" % (missing_desc, response))

    def __init__(self,
            db,
            ids=None,
            id_range=None,
            prefix=None,
            verify=True,
            email=None):
        """
        Initializes a GenBank data resource manager class, and optionally
        populates it with data from GenBank.

            - ``db``: database (e.g. "nucleotide" or "protein")
            - ``ids``: sequence of GenBank accession identifiers or GI numbers.
            - ``id_range``: tuple of integers indicating first and last (inclusive) range of identifiers to download.
            - ``prefix``: string to be prepended to values in - ``ids``: or generated from - ``id_range``:
            - ``verify``: check for one-to-one correspondence between requested id's and downloaded records

        """
        self.db = db
        self.email = email
        self._recs = []
        self._accession_recs = {}
        self._version_recs = {}
        self._gi_recs = {}
        if ids is not None:
            self.acquire(ids=ids,
                    prefix=prefix,
                    verify=verify)
        if id_range is not None:
            self.acquire_range(
                    first=id_range[0],
                    last=id_range[1],
                    prefix=prefix,
                    verify=verify)

    def __len__(self):
        return len(self._recs)

    def __getitem__(self, key):
        return self._recs[key]

    def __setitem__(self, key, value):
        raise TypeError("%s elements cannot be reassigned" % self.__class__)

    def __delitem__(self, key):
        raise TypeError("%s elements cannot be deleted" % self.__class__)

    def __iter__(self):
        return iter(self._recs)

    def __reversed__(self):
        return reversed(self._recs)

    def __contains__(self, item):
        return item in self._recs

    def add(self, rec):
        if rec.primary_accession in self._accession_recs:
            return None
        if rec.accession_version in self._version_recs:
            return None
        if rec.gi in self._gi_recs:
            return None
        if rec in self._recs:
            return None
        if rec.db is None:
            rec.db = self.db
        self._recs.append(rec)
        self._accession_recs[rec.primary_accession] = rec
        self._version_recs[rec.accession_version] = rec
        self._gi_recs[rec.gi] = rec

    def update(self, recs):
        for rec in recs:
            self.add(rec)

    def read_xml_string(self, xml_string):
        gb_recs = GenBankResourceStore.parse_xml(string=xml_string)
        self.update(gb_recs)

    def read_xml_stream(self, xml_stream):
        gb_recs = GenBankResourceStore.parse_xml(stream=xml_stream)
        self.update(gb_recs)

    def acquire(self, ids, prefix=None, verify=True):
        """
        Adds more data from GenBank.

            - ``ids``: sequence of GenBank accession identifiers or GI numbers.
            - ``prefix``: string to be prepended to values in - ``ids`` or generated from - ``id_range``.
            - ``verify``: check for one-to-one correspondence between requested id's and downloaded records

        """
        ids = GenBankResourceStore.prepare_ids(ids=ids, prefix=prefix)
        xml_string = GenBankResourceStore.fetch_xml(db=self.db,
                ids=ids,
                prefix=prefix,
                email=self.email,
                as_stream=False)
        gb_recs = GenBankResourceStore.parse_xml(string=xml_string)
        accession_recs = {}
        accession_version_recs = {}
        gi_recs = {}
        for gb_rec in gb_recs:
            accession_recs[gb_rec.primary_accession] = gb_rec
            accession_version_recs[gb_rec.accession_version] = gb_rec
            gi_recs[gb_rec.gi] = gb_rec
        added = []
        missing = []
        for idx, gbid in enumerate(ids):
            sgbid = str(gbid)
            gb_rec = None
            if sgbid in accession_version_recs:
                gb_rec = accession_version_recs[sgbid]
            elif sgbid in accession_recs:
                gb_rec = accession_recs[sgbid]
            elif sgbid in gi_recs:
                gb_rec = gi_recs[sgbid]
            elif verify:
                missing.append(sgbid)
            if gb_rec is not None:
                gb_rec.db = self.db
                gb_rec.request_key = sgbid
                added.append(gb_rec)
        if len(added) == 0 and missing:
            raise GenBankResourceStore.AccessionFetchError(missing=missing, response=xml_string)
        elif missing:
            raise GenBankResourceStore.AccessionFetchError(missing=missing, response=None)
        self.update(added)
        return added

    def acquire_range(self,
            first,
            last,
            prefix=None,
            verify=True):
        """
        Adds more data from GenBank.

            - ``first``: integer specifying the start (inclusive) of the range of id's to download.
            - ``last``: integer specifying the start (inclusive) of the range of id's to download.
            - ``prefix``: string to be prepended to values given range.
            - ``verify``: check for one-to-one correspondence between requested id's and downloaded records

        """
        ids = range(first, last+1)
        return self.acquire(
                ids=ids,
                prefix=prefix,
                verify=verify)

    def generate_char_matrix(self,
            label_components=None,
            label_component_separator=" ",
            taxon_namespace=None,
            gb_to_taxon_fn=None,
            add_full_annotation_to_taxa=False,
            add_ref_annotation_to_taxa=False,
            add_full_annotation_to_seqs=False,
            add_ref_annotation_to_seqs=False,
            set_taxon_attr=None,
            set_seq_attr=None,
            matrix_label=None):
        """
        Generates a CharacterMatrix of current sequences.

            - ``label_components``: list of strings giving names of GenBankAccessionRecord attributes to be used to compose label.
            - ``label_component_separator``: a string used to separate label components.
            - ``taxon_namespace``: TaxonNamespace object to be used as the ``taxon_namespace`` of the resulting CharacterMatrix.
            - ``gb_to_taxon_fn``: Function to be used to assign a Taxon object to sequence. Should take a GenBankAccessionRecord object as an argument and return a Taxon object.
            - ``add_full_annotation_to_taxa``: If True, add link to record as metadata annotation to Taxon objects.
            - ``add_ref_annotation_to_taxa``: If True, add full GenBank record as metadata annotations to Taxon objects.
            - ``add_full_annotation_to_seqs``: If True, add link to record as metadata annotation to sequence (CharacterDataSequence) objects.
            - ``add_ref_annotation_to_seqs``: If True, add full GenBank record as metadata annotations to sequence (CharacterDataSequence) objects.
            - ``set_taxon_attr``: Name of attribute (string) to create on Taxon objects pointing to GenBank record object (GenBankAccessionRecord).
            - ``set_seq_attr``: Name of attribute (string) to create on sequence objects pointing to GenBank record object (GenBankAccessionRecord).
            - ``matrix_label``: Label of character matrix.

        """
        if gb_to_taxon_fn is not None and taxon_namespace is None:
            raise TypeError("Cannot specify 'gb_to_taxon_fn' without 'taxon_namespace'")
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        data_str = []
        char_matrix = self.char_matrix_type(label=matrix_label, taxon_namespace=taxon_namespace)
        for gb_idx, gb_rec in enumerate(self._recs):
            taxon = None
            # if gb_rec.request_key in id_to_taxon_map:
            #     taxon = id_to_taxon_map[gb_rec.request_key]
            if gb_to_taxon_fn is not None:
                taxon = gb_to_taxon_fn(gb_rec)
            else:
                if label_components is not None:
                    label = gb_rec.compose_taxon_label(
                            components=label_components,
                            separator=label_component_separator)
                else:
                    label = gb_rec.request_key
                assert label is not None
                assert str(label) != "None"
                taxon = taxon_namespace.require_taxon(label=label)
            assert taxon is not None
            if add_ref_annotation_to_taxa:
                taxon.annotations.add(gb_rec.as_reference_annotation())
            if add_full_annotation_to_taxa:
                taxon.annotations.add(gb_rec.as_annotation())
            if set_taxon_attr is not None:
                setattr(taxon, set_taxon_attr, gb_rec)
            curr_vec = char_matrix.new_sequence(taxon=taxon)
            char_matrix[taxon] = curr_vec
            if add_ref_annotation_to_seqs:
                curr_vec.annotations.add(gb_rec.as_reference_annotation())
            if add_full_annotation_to_seqs:
                curr_vec.annotations.add(gb_rec.as_annotation())
            if set_seq_attr is not None:
                setattr(curr_vec, set_seq_attr, gb_rec)
            seq_text = gb_rec.sequence_text.upper()
            for col_ind, c in enumerate(seq_text):
                c = c.strip()
                if not c:
                    continue
                try:
                    state = char_matrix.default_state_alphabet[c]
                except KeyError:
                    raise ValueError('Accession %d of %d (%s, GI %s, acquired using key: %s): Unrecognized sequence symbol "%s" in position %d: %s'
                            % (gb_idx+1, len(self._recs), gb_rec.primary_accession, gb_rec.gi, gb_rec.request_key, c, col_ind+1, seq_text))
                curr_vec.append(state)
        return char_matrix

class GenBankNucleotide(GenBankResourceStore):

    def __init__(self,
            ids=None,
            id_range=None,
            prefix=None,
            verify=True,
            char_matrix_type=None,
            email=None):
        GenBankResourceStore.__init__(self,
                db="nucleotide",
                ids=ids,
                id_range=id_range,
                prefix=prefix,
                verify=verify,
                email=email)
        self.char_matrix_type = char_matrix_type

class GenBankDna(GenBankNucleotide):

    def __init__(self,
            ids=None,
            id_range=None,
            prefix=None,
            verify=True,
            email=None):
        GenBankNucleotide.__init__(self,
                ids=ids,
                id_range=id_range,
                prefix=prefix,
                verify=verify,
                char_matrix_type=dendropy.DnaCharacterMatrix,
                email=email)

class GenBankRna(GenBankNucleotide):

    def __init__(self,
            ids=None,
            id_range=None,
            prefix=None,
            verify=True,
            email=None):
        GenBankNucleotide.__init__(self,
                ids=ids,
                id_range=id_range,
                prefix=prefix,
                verify=verify,
                char_matrix_type=dendropy.RnaCharacterMatrix,
                email=email)

class GenBankProtein(GenBankResourceStore):

    def __init__(self,
            ids=None,
            id_range=None,
            prefix=None,
            verify=True,
            email=None):
        GenBankResourceStore.__init__(self,
                db="protein",
                ids=ids,
                id_range=id_range,
                prefix=prefix,
                verify=verify,
                email=email)
        self.char_matrix_type = dendropy.ProteinCharacterMatrix

##############################################################################
## GenBank Data Parsing to Python Objects

class GenBankAccessionReference(object):
    """
    A GenBank reference record.
    """
    def __init__(self, xml=None):
        self.number = None
        self.position = None
        self.authors = []
        self.consrtm = None
        self.title = None
        self.journal = None
        self.medline_id = None
        self.pubmed_id = None
        self.remark = None
        self.db_ids = container.OrderedCaselessDict()
        if xml is not None:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.number = xml.findtext("INSDReference_reference")
        self.position = xml.findtext("INSDReference_position")
        for authors in xml.findall("INSDReference_authors"):
            for author_xml in authors.findall("INSDReference_author"):
                author = author_xml.findtext("INSDAuthor")
                self.authors.append(author)
        self.title = xml.findtext("INSDReference_title")
        self.journal = xml.findtext("INSDReference_journal")
        for xrefs in xml.findall("INSDReference_xref"):
            for xref_xml in xrefs.findall("INSDXref"):
                dbname = xref_xml.findtext("INSDXref_dbname")
                dbid = xref_xml.findtext("INSDXref_id")
                self.db_ids[dbname] = dbid
            self.pubmed_id = xml.findtext("INSDReference_pubmed")
            self.medline_id = xml.findtext("INSDReference_medline")

    def as_annotation(self):
        top = dendropy.Annotation(
                name="INSDReference_reference",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for item in [
                ("number", "INSDReference_reference"),
                ("position", "INSDReference_position"),
                ("title", "INSDReference_title"),
                ("journal", "INSDReference_journal"),
                ("pubmed_id", "INSDReference_pubmed"),
                ("medline_id", "INSDReference_medline"),
                ]:
            value = getattr(self, item[0])
            if value is None:
                continue
            sub = dendropy.Annotation(
                name=item[1],
                value=value,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=False,
                is_hidden=False,
                label=None,
                oid=None)
            top.annotations.add(sub)
        return top


class GenBankAccessionReferences(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def findall(self, key):
        results = []
        for a in self:
            if a.key == key:
                results.append(a)
        results = GenBankAccessionFeatures(results)
        return results

    def find(self, key, default=None):
        for a in self:
            if a.key == key:
                return a
        return default

    def get_value(self, key, default=None):
        for a in self:
            if a.key == key:
                return a.value
        return default

    def as_annotation(self):
        top = dendropy.Annotation(
                name="INSDSeq_references",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for reference in self:
            top.annotations.add(reference.as_annotation())
        return top

class GenBankAccessionInterval(object):
    """
    A GenBank Interval record.
    """
    def __init__(self, xml=None):
        self.begin = None
        self.end = None
        self.accession = None
        if xml is not None:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.begin = xml.findtext("INSDInterval_from")
        self.end = xml.findtext("INSDInterval_to")
        self.accession = xml.findtext("INSDInterval_accession")

class GenBankAccessionQualifier(object):
    """
    A GenBank Qualifier record.
    """
    def __init__(self, xml=None):
        self.name = None
        self.value = None
        if xml is not None:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.name = xml.findtext("INSDQualifier_name")
        self.value = xml.findtext("INSDQualifier_value")

    def as_annotation(self):
        annote = dendropy.Annotation(
            name=self.name,
            value=self.value,
            datatype_hint=None,
            name_prefix=GENBANK_ANNOTATION_PREFIX,
            namespace=GENBANK_ANNOTATION_NAMESPACE,
            name_is_prefixed=False,
            is_attribute=False,
            annotate_as_reference=False,
            is_hidden=False,
            label=None,
            oid=None)
        return annote

class GenBankAccessionQualifiers(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def findall(self, name):
        results = []
        for a in self:
            if a.name == name:
                results.append(a)
        results = GenBankAccessionQualifiers(results)
        return results

    def find(self, name, default=None):
        for a in self:
            if a.name == name:
                return a
        return default

    def get_value(self, name, default=None):
        for a in self:
            if a.name == name:
                return a.value
        return default

    def as_annotation(self):
        top = dendropy.Annotation(
                name="INSDFeature_quals",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for item in self:
            top.annotations.add(item.as_annotation())
        return top

class GenBankAccessionFeature(object):
    """
    A GenBank Feature record.
    """
    def __init__(self, xml=None):
        self.key = None
        self.location = None
        self.qualifiers = GenBankAccessionQualifiers()
        self.intervals = []
        if xml is not None:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.key = xml.findtext("INSDFeature_key")
        self.location = xml.findtext("INSDFeature_location")
        for intervals in xml.findall("INSDFeature_intervals"):
            for interval_xml in xml.findall("INSDInterval"):
                interval = GenBankAccessionInterval(interval_xml)
                self.intervals.append(interval)
        for qualifiers in xml.findall("INSDFeature_quals"):
            for qualifier_xml in qualifiers.findall("INSDQualifier"):
                qualifier = GenBankAccessionQualifier(qualifier_xml)
                self.qualifiers.append(qualifier)

    def as_annotation(self):
        top = dendropy.Annotation(
                name="INSDSeq_feature",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for item in [
                ("key", "INSDFeature_key"),
                ("location", "INSDFeature_location"),
                ]:
            value = getattr(self, item[0])
            if not value:
                continue
            sub = dendropy.Annotation(
                name=item[1],
                value=value,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=False,
                is_hidden=False,
                label=None,
                oid=None)
            top.annotations.add(sub)
        if self.intervals:
            intervals_annote = dendropy.Annotation(
                name="INSDSeq_intervals",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
            top.annotations.add(intervals_annote)
            for interval in self.intervals:
                interval_annote = dendropy.Annotation(
                    name="INSDInterval",
                    value=None,
                    datatype_hint=None,
                    name_prefix=GENBANK_ANNOTATION_PREFIX,
                    namespace=GENBANK_ANNOTATION_NAMESPACE,
                    name_is_prefixed=False,
                    is_attribute=False,
                    annotate_as_reference=True,
                    is_hidden=False,
                    label=None,
                    oid=None)
                intervals_annote.annotations.add(interval_annote)
                for item in [
                        ("begin", "INSDInterval_from"),
                        ("end", "INSDInterval_to"),
                        ("accession", "INSDInterval_accession"),
                        ]:
                    value = getattr(self, item[0])
                    if not value:
                        continue
                    sub = dendropy.Annotation(
                        name=item[1],
                        value=value,
                        datatype_hint=None,
                        name_prefix=GENBANK_ANNOTATION_PREFIX,
                        namespace=GENBANK_ANNOTATION_NAMESPACE,
                        name_is_prefixed=False,
                        is_attribute=False,
                        annotate_as_reference=False,
                        is_hidden=False,
                        label=None,
                        oid=None)
                    interval_annote.annotations.add(interval_annote)
        if self.qualifiers:
            top.annotations.add(self.qualifiers.as_annotation())
        return top

class GenBankAccessionFeatures(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def findall(self, key):
        results = []
        for a in self:
            if a.key == key:
                results.append(a)
        results = GenBankAccessionFeatures(results)
        return results

    def find(self, key, default=None):
        for a in self:
            if a.key == key:
                return a
        return default

    def get_value(self, key, default=None):
        for a in self:
            if a.key == key:
                return a.value
        return default

    def as_annotation(self):
        top = dendropy.Annotation(
                name="INSDSeq_feature-table",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for feature in self:
            top.annotations.add(feature.as_annotation())
        return top

class GenBankAccessionOtherSeqIds(dict):

    def __init__(self, *args):
        dict.__init__(self, *args)

    def as_annotation(self):
        top = dendropy.Annotation(
                name="otherSeqIds",
                value=None,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        for key, value in self.items():
            sub = dendropy.Annotation(
                name=key,
                value=value,
                datatype_hint=None,
                name_prefix=GENBANK_ANNOTATION_PREFIX,
                namespace=GENBANK_ANNOTATION_NAMESPACE,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=False,
                is_hidden=False,
                label=None,
                oid=None)
            top.annotations.add(sub)
        return top

class GenBankAccessionRecord(object):
    """
    A GenBank record.
    """

    def __init__(self, db=None, xml=None):
        self.db = db
        self._request_key = None
        self._defline = None
        self.locus = None
        self.length = None
        self.moltype = None
        self.topology = None
        self.strandedness = None
        self.division = None
        self.update_date = None
        self.create_date = None
        self.definition = None
        self.primary_accession = None
        self.accession_version = None
        self.other_seq_ids = GenBankAccessionOtherSeqIds()
        self.source = None
        self.organism = None
        self.taxonomy = None
        self.references = GenBankAccessionReferences()
        self.feature_table = GenBankAccessionFeatures()
        self.sequence_text = None
        if xml is not None:
            self.parse_xml(xml)

    def __str__(self):
        return self.defline

    def _get_defline(self):
        if self._defline is None:
            self._defline = self.compose_fasta_defline()
        return self._defline
    defline = property(_get_defline)

    def _get_request_key(self):
        if self._request_key is None:
            return self.gi
        else:
            return self._request_key
    def _set_request_key(self, rid):
        self._request_key = rid
    request_key = property(_get_request_key, _set_request_key)

    def _has_request_key(self):
        return self._request_key is not None
    has_request_key = property(_has_request_key)

    def _get_accession(self):
        return self.primary_accession
    accession = property(_get_accession)

    def _get_gi(self):
        return self.other_seq_ids.get("gi", None)
    gi = property(_get_gi)

    def _get_gb(self):
        return self.other_seq_ids.get("gb", None)
    gb = property(_get_gb)

    def parse_xml(self, xml):
        self.locus = xml.findtext("INSDSeq_locus")
        self.length = xml.findtext("INSDSeq_length")
        self.moltype = xml.findtext("INSDSeq_moltype")
        self.topology = xml.findtext("INSDSeq_topology")
        self.strandedness = xml.findtext("INSDSeq_strandedness")
        self.division = xml.findtext("INSDSeq_division")
        self.update_date = xml.findtext("INSDSeq_update-date")
        self.create_date = xml.findtext("INSDSeq_create-date")
        self.definition = xml.findtext("INSDSeq_definition")
        self.primary_accession = xml.findtext("INSDSeq_primary-accession")
        self.accession_version = xml.findtext("INSDSeq_accession-version")
        for other_seqids in xml.findall("INSDSeq_other-seqids"):
            for other_seqid in other_seqids.findall("INSDSeqid"):
                seqid = other_seqid.text
                parts = seqid.split("|")
                if len(parts) == 1:
                    self.other_seq_ids[""] = seqid
                else:
                    self.other_seq_ids[parts[0]] = parts[1]
        self.source = xml.findtext("INSDSeq_source")
        self.organism = xml.findtext("INSDSeq_organism")
        self.taxonomy = xml.findtext("INSDSeq_taxonomy")
        for references in xml.findall("INSDSeq_references"):
            for reference_xml in references.findall("INSDReference"):
                reference = GenBankAccessionReference(reference_xml)
                self.references.append(reference)
        for features in xml.findall("INSDSeq_feature-table"):
            for feature_xml in features.findall("INSDFeature"):
                feature = GenBankAccessionFeature(feature_xml)
                self.feature_table.append(feature)
        self.sequence_text = xml.findtext("INSDSeq_sequence")

    def compose_fasta_defline(self):
        return "gi|%s|gb|%s| %s" % (
                self.gi,
                self.accession_version,
                self.definition)

    def compose_taxon_label(self, components=None, separator=" "):
        label = []
        if components is None:
            components = ["accession", "organism"]
        for component in components:
            component = str(getattr(self, component))
            if component is not None:
                if separator is not None and " " in component and separator != " ":
                    component = component.replace(" ", separator)
                label.append(component)
        if separator is None:
            separator = ""
        return separator.join(label)

    def as_fasta(self,
            generate_new_label=False,
            label_components=None,
            label_component_separator=" "):
        if generate_new_label:
            label = self.compose_taxon_label(components=label_components,
                    separator=label_component_separator)
        else:
            label = self.compose_fasta_defline()
        sequence_text = self.sequence_text.upper()
        return ">%s\n%s" % (label, sequence_text)

    def _get_uri(self):
        uri = ["http://www.ncbi.nlm.nih.gov"]
        uri.append(str(self.db))
        uri.append(str(self.request_key))
        return "/".join(uri)
    uri = property(_get_uri)

    def as_reference_annotation(self):
        annote = dendropy.Annotation(
                name="source",
                value=self.uri,
                datatype_hint=None,
                name_prefix="dcterms",
                namespace="http://purl.org/dc/terms/",
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                is_hidden=False,
                label=None,
                oid=None)
        return annote

    def as_annotation(self):
        top = self.as_reference_annotation()
        # top = dendropy.Annotation(
        #         name="source",
        #         value=None,
        #         datatype_hint=None,
        #         name_prefix="dcterms",
        #         namespace="http://purl.org/dc/terms/",
        #         name_is_prefixed=False,
        #         is_attribute=False,
        #         annotate_as_reference=True,
        #         is_hidden=False,
        #         label=None,
        #         oid=None)
        # for item in [
        #         ("db", "db"),
        #         ("request_key", "requestKey")
        #         ]:
        #     value = getattr(self, item[0])
        #     if value is None:
        #         continue
        #     a = dendropy.Annotation(
        #         name=item[1],
        #         value=value,
        #         datatype_hint=None,
        #         name_prefix="dendropy",
        #         namespace="http://packages.python.org/DendroPy/",
        #         name_is_prefixed=False,
        #         is_attribute=False,
        #         annotate_as_reference=False,
        #         is_hidden=False,
        #         label=None,
        #         oid=None)
        #     top.annotations.add(a)
        for item in [
                ("locus", "INSDSeq_locus"),
                ("length", "INSDSeq_length"),
                ("moltype", "INSDSeq_moltype"),
                ("topology", "INSDSeq_topology"),
                ("strandedness", "INSDSeq_strandedness"),
                ("division", "INSDSeq_division"),
                ("update_date", "INSDSeq_update-date"),
                ("create_date", "INSDSeq_create-date"),
                ("definition", "INSDSeq_definition"),
                ("primary_accession", "INSDSeq_primary-accesison"),
                ("accession_version", "INSDSeq_accession-version"),
                ("other_seq_ids",  "INSDSeq_hter-seqids"),
                ("source", "INSDSeq_source"),
                ("organism", "INSDSeq_organism"),
                ("taxonomy", "INSDSeq_taxonomy"),
                ("references", "INSDSeq_references"),
                ("feature_table", "INSDSeq_feature-table"),
                ]:
            value = getattr(self, item[0])
            if hasattr(value, "as_annotation") and value:
                a = value.as_annotation()
            elif value is not None:
                a = dendropy.Annotation(
                    name=item[1],
                    value=value,
                    datatype_hint=None,
                    name_prefix=GENBANK_ANNOTATION_PREFIX,
                    namespace=GENBANK_ANNOTATION_NAMESPACE,
                    name_is_prefixed=False,
                    is_attribute=False,
                    annotate_as_reference=False,
                    is_hidden=False,
                    label=None,
                    oid=None)
            top.annotations.add(a)
        return top

