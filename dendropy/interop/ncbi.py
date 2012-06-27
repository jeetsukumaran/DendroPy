#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Wrappers for interacting with NCBI databases.
"""

import re
import urllib
import sys
import dendropy
from dendropy.dataio.xmlparser import ElementTree
from dendropy.utility.error import DataParseError
from dendropy.utility import containers
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

# GB_FASTA_DEFLINE_PATTERN = re.compile(r'^gi\|(\d+)\|(\w+)\|([\w\d]+).(\d+)\|(.*)$')

# def parse_accession_number_and_gi_from_gb(gb_str):
#     accession = re.search(r"ACCESSION\s+(\S+)$", gb_str, flags=re.MULTILINE)
#     if accession is None:
#         raise ValueError("Failed to parse accession number")
#     accession = accession.groups()[0].strip()
#     gi = re.search(r'^VERSION\s+\S+\s+GI:([0-9]+)$', gb_str, flags=re.MULTILINE)
#     if gi is None:
#         raise ValueError("Failed to parse GI number")
#     gi = gi.groups()[0].strip()
#     return accession, gi

# def parse_ncbi_curation_info_from_defline(gb_defline):
#     m = GB_FASTA_DEFLINE_PATTERN.match(gb_defline)
#     if m is not None:
#         return m.groups()[0], m.groups()[2], m.groups()[2] + '.' + m.groups()[3]
#     else:
#         return None

# def compose_taxon_label_from_gb_defline(gb_defline,
#         num_desc_components=3,
#         separator='_',
#         gbnum_in_front=True,
#         exclude_gbnum=False):
#     """
#     If `gb_defline` matches a GenBank FASTA-format defline structure, then this returns a
#     label:

#         <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

#     So, for example, given the following FASTA label:

#         gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

#     the corresponding taxon 3-component (default) label will be:

#         EU105975_Homo_sapiens_Ache

#     If `gb_defline` does *not* match a GenBank FASTA-format defline structure, then the string
#     is returned unchanged.
#     """
#     m = GB_FASTA_DEFLINE_PATTERN.match(gb_defline)
#     if m is not None:
#         groups = m.groups()
#         desc_parts = [s.strip() for s in groups[-1].split() if s]
#         if exclude_gbnum:
#             label_parts = desc_parts[:num_desc_components]
#         elif gbnum_in_front:
#             label_parts = [groups[2]] + desc_parts[:num_desc_components]
#         else:
#             label_parts = desc_parts[:num_desc_components] + [groups[2]]
#         return separator.join(label_parts)
#     else:
#         return gb_defline

# def relabel_taxa_from_defline(taxon_set,
#         num_desc_components=3,
#         separator='_',
#         gbnum_in_front=True,
#         exclude_gbnum=False):
#     """
#     Examines the labels of each `Taxon` object in `taxon_set`, and if
#     conforming to a GenBank pattern, translates the labels to a standard
#     format:

#         <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

#     So, for example, given the following FASTA label:

#         gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

#     the corresponding taxon 3-component (default) label will be:

#         EU105975_Homo_sapiens_Ache

#     """
#     for taxon in taxon_set:
#         taxon.label = compose_taxon_label_from_gb_defline(
#                 gb_defline=taxon.label,
#                 num_desc_components=num_desc_components,
#                 separator=separator,
#                 gbnum_in_front=gbnum_in_front,
#                 exclude_gbnum=exclude_gbnum)
#     return taxon_set

class GenBankReference(object):
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
        self.db_ids = containers.OrderedCaselessDict()
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

class GenBankInterval(object):
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

class GenBankQualifier(object):
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

class GenBankQualifiers(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def findall(self, name):
        results = []
        for a in self:
            if a.name == name:
                results.append(a)
        results = GenBankQualifiers(results)
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

class GenBankFeature(object):
    """
    A GenBank Feature record.
    """
    def __init__(self, xml=None):
        self.key = None
        self.location = None
        self.qualifiers = GenBankQualifiers()
        self.intervals = []
        if xml is not None:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.key = xml.findtext("INSDFeature_key")
        self.location = xml.findtext("INSDFeature_location")
        for intervals in xml.findall("INSDFeature_intervals"):
            for interval_xml in xml.findall("INSDInterval"):
                interval = GenBankInterval(interval_xml)
                self.intervals.append(interval)
        for qualifiers in xml.findall("INSDFeature_quals"):
            for qualifier_xml in qualifiers.findall("INSDQualifier"):
                qualifier = GenBankQualifier(qualifier_xml)
                self.qualifiers.append(qualifier)

class GenBankFeatures(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def findall(self, key):
        results = []
        for a in self:
            if a.key == key:
                results.append(a)
        results = GenBankFeatures(results)
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

class GenBankRecord(object):
    """
    A GenBank record.
    """

    def parse_from_stream(stream):
        tree = ElementTree.parse(stream)
        root = tree.getroot()
        gb_recs = []
        for seq_set in root.iter("INSDSet"):
            for seq in seq_set.iter("INSDSeq"):
                gb_rec = GenBankRecord(xml=seq)
                gb_recs.append(gb_rec)
        return gb_recs
    parse_from_stream = staticmethod(parse_from_stream)

    def __init__(self, xml=None):
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
        self.other_seq_ids = {}
        self.source = None
        self.organism = None
        self.taxonomy = None
        self.references = []
        self.features = GenBankFeatures()
        self.sequence_text = None
        if xml is not None:
            self.parse_xml(xml)

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
                reference = GenBankReference(reference_xml)
                self.references.append(reference)
        for features in xml.findall("INSDSeq_feature-table"):
            for feature_xml in features.findall("INSDFeature"):
                feature = GenBankFeature(feature_xml)
                self.features.append(feature)
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


class Entrez(object):
    """
    Wraps up all interactions with Entrez.
    Example usage::

        >>> from dendropy.interop import ncbi
        >>> e = ncbi.Entrez(generate_labels=True,
        ... label_id_in_front=False,
        ... sort_taxa_by_label=True)
        >>> d1 = e.fetch_nucleotide_accessions(['EU105474', 'EU105476'])
        >>> d2 = e.fetch_nucleotide_accession_range(105474, 106045, prefix="EU")

    """

    BASE_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    DATABASES = [
        'pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest', 'structure',
        'genome', 'biosystems', 'books', 'cancerchromosomes', 'cdd', 'gap', 'dbvar',
        'domains', 'epigenomics', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
        'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pepdome', 'pmc',
        'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound', 'pcsubstance', 'seqannot',
        'snp', 'sra', 'taxonomy', 'toolkit', 'toolkitall', 'unigene', 'unists', 'linkoutpubmed',
        'linkoutseq', 'linkoutother',
        ]

    class AccessionFetchError(Exception):

        def __init__(self, accession_ids):
            Exception.__init__(self, "Failed to retrieve accessions: %s" % (", ".join([str(s) for s in accession_ids])))

    def __init__(self):
        """
        Instantiates a broker that queries NCBI and returns data.
        """
        pass

    def efetch(self, db, ids, rettype, retmode="text"):
        """
        Raw fetch. Returns file-like object opened for reading on string
        returned by query.
        """
        if isinstance(ids, str):
            id_list = ids
        else:
            id_list = ",".join([str(i) for i in set(ids)])
        params = {'db': db,
                'id': id_list,
                'rettype': rettype,
                'retmode': retmode}
        query_url = Entrez.BASE_URL + "/efetch.fcgi?" + urllib.urlencode(params)
        query = urllib.urlopen(query_url)
        return query

#     def fetch_gbrecs_as_plaintext_dict(self, db, ids, verify=True):
#         db_name = "nucleotide"
#         gb_recs_str = self.fetch(db=db, ids=ids, rettype="gb")
#         gb_recs_str_list = re.split(r"^//$", gb_recs_str, flags=re.MULTILINE)
#         gb_recs_str_list = [gb_rec for gb_rec in gb_recs_str_list
#                     if gb_rec.replace("\n", "")]
#         accession_recs = {}
#         gi_recs = {}
#         for gb_str in gb_recs_str_list:
#             if not gb_str:
#                 continue
#             accession, gi = parse_accession_number_and_gi_from_gb(gb_str)
#             accession_recs[accession] = gb_str
#             gi_recs[str(gi)] = gb_str
#         result = containers.OrderedCaselessDict()
#         for gbid in ids:
#             sgbid = str(gbid)
#             if sgbid in accession_recs:
#                 result[gbid] = accession_recs[sgbid]
#             elif sgbid in gi_recs:
#                 result[gbid] = gi_recs[sgbid]
#             elif verify:
#                 raise Entrez.AccessionFetchError(sgbid)
#         return result

    def fetch_genbank_records(self, db, ids, verify=True):
        results_stream = self.efetch(db=db,
                ids=ids,
                rettype='gbc',
                retmode='xml')
        gb_recs = GenBankRecord.parse_from_stream(results_stream)
        accession_recs = {}
        accession_version_recs = {}
        gi_recs = {}
        for gb_rec in gb_recs:
            accession_recs[gb_rec.accession] = gb_rec
            accession_version_recs[gb_rec.accession_version] = gb_rec
            gi_recs[gb_rec.gi] = gb_rec
        result = []
        missing = []
        for idx, gbid in enumerate(ids):
            sgbid = str(gbid)
            if sgbid in accession_version_recs:
                result.append( accession_version_recs[sgbid] )
            elif sgbid in accession_recs:
                result.append( accession_recs[sgbid] )
            elif sgbid in gi_recs:
                result.append( gi_recs[sgbid] )
            elif verify:
                missing.append(sgbid)
        if missing:
            raise Entrez.AccessionFetchError(missing)
        return result

    def fetch_nucleotide_accessions(self,
            ids,
            prefix=None,
            verify=True,
            matrix_type=dendropy.DnaCharacterMatrix,
            generate_new_labels=False,
            label_components=None,
            label_component_separator=" ",
            sort_taxa_by_label=False,
            **kwargs):
        """
        Returns a DnaCharacterMatrix object (or some other type, if specified
        by ``matrix_type`` argument) populated with sequences from the Entrez
        nucleotide database with accession numbers given by `ids` (a list of
        accession numbers). If `prefix` is given, it is pre-pended to all values
        given in the id list. Any other keyword arguments given are passed to
        the constructor of ``DnaCharacterMatrix``.
        By default, the sequence labels are as given by the GenBank defline.
        If ``generate_new_labels`` is True then the sequence labels will be
        generated based on the names of the GenBank attributes passed in
        label components. This defaults to ``["accession", "organism"]``.
        """
        if prefix is not None:
            ids = ["%s%s" % (prefix,i) for i in ids]
        # results_stream = self.efetch(db='nucleotide',
        #         ids=ids,
        #         rettype='gbc',
        #         retmode='xml')
        gb_recs = self.fetch_genbank_records(db="nucleotide",
                ids=ids,
                verify=verify)
        data_str = []
        for gb_rec in gb_recs:
            fasta = gb_rec.as_fasta(generate_new_label=generate_new_labels,
                    label_components=label_components,
                    label_component_separator=label_component_separator)
            data_str.append(fasta)
        data_str = "\n".join(data_str)
        d = matrix_type.get_from_string(data_str,
                'fasta',
                **kwargs)
        if sort_taxa_by_label:
            d.taxon_set.sort(key=lambda x: x.label)
        return d

    def fetch_nucleotide_accession_range(self,
            first,
            last,
            prefix=None,
            verify=True,
            matrix_type=dendropy.DnaCharacterMatrix,
            **kwargs):
        """
        Returns a DnaCharacterMatrix object (or some other type, if specified
        by the ``matrix_type`` argument) populated with sequences from the
        Entrez nucleotide database with accession numbers between ``start``
        and, up to and *including* ``end``. If `prefix` is given, then it is
        pre-pended to the ids. Any other keyword arguments given are passed to
        thee constructor of ``DnaCharacterMatrix``.
        """
        ids = range(first, last+1)
        return self.fetch_nucleotide_accessions(ids=ids, prefix=prefix, verify=verify, matrix_type=matrix_type, **kwargs)

