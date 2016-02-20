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
Wrappers for interacting with NCBI databases.
*** DEPRECATED: use dendropy.interop.genbank.GenBankDna,
dendropy.interop.genbank.GenBankRna, or
dendropy.interop.genbank.GenBankProtein instead ***
"""

import warnings
from dendropy.utility import messaging
from dendropy.utility import textprocessing
from dendropy.utility import urlio
_LOG = messaging.get_logger(__name__)

import sys
import dendropy
import re
from dendropy.utility.error import DataParseError

GB_FASTA_DEFLINE_PATTERN = re.compile(r'^gi\|(\d+)\|(\w+)\|([\w\d]+).(\d+)\|(.*)$')

def parse_accession_number_and_gi_from_gb(gb_str):
    accession = re.search(r"ACCESSION\s+(\S+)$", gb_str, flags=re.MULTILINE)
    if accession is None:
        raise ValueError("Failed to parse accession number")
    accession = accession.groups()[0].strip()
    gi = re.search(r'^VERSION\s+\S+\s+GI:([0-9]+)$', gb_str, flags=re.MULTILINE)
    if gi is None:
        raise ValueError("Failed to parse GI number")
    gi = gi.groups()[0].strip()
    return accession, gi

def parse_ncbi_curation_info_from_defline(gb_defline):
    m = GB_FASTA_DEFLINE_PATTERN.match(gb_defline)
    if m is not None:
        return m.groups()[0], m.groups()[2], m.groups()[2] + '.' + m.groups()[3]
    else:
        return None

def compose_taxon_label_from_gb_defline(gb_defline,
        num_desc_components=3,
        separator='_',
        gbnum_in_front=True,
        exclude_gbnum=False):
    """
    If ``gb_defline`` matches a GenBank FASTA-format defline structure, then this returns a
    label:

        <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

    So, for example, given the following FASTA label:

        gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

    the corresponding taxon 3-component (default) label will be:

        EU105975_Homo_sapiens_Ache

    If ``gb_defline`` does *not* match a GenBank FASTA-format defline structure, then the string
    is returned unchanged.
    """
    m = GB_FASTA_DEFLINE_PATTERN.match(gb_defline)
    if m is not None:
        groups = m.groups()
        desc_parts = [s.strip() for s in groups[-1].split() if s]
        if exclude_gbnum:
            label_parts = desc_parts[:num_desc_components]
        elif gbnum_in_front:
            label_parts = [groups[2]] + desc_parts[:num_desc_components]
        else:
            label_parts = desc_parts[:num_desc_components] + [groups[2]]
        return separator.join(label_parts)
    else:
        return gb_defline

def relabel_taxa_from_defline(taxon_set,
        num_desc_components=3,
        separator='_',
        gbnum_in_front=True,
        exclude_gbnum=False):
    """
    Examines the labels of each |Taxon| object in ``taxon_set``, and if
    conforming to a GenBank pattern, translates the labels to a standard
    format:

        <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

    So, for example, given the following FASTA label:

        gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

    the corresponding taxon 3-component (default) label will be:

        EU105975_Homo_sapiens_Ache

    """
    for taxon in taxon_set:
        taxon.label = compose_taxon_label_from_gb_defline(
                gb_defline=taxon.label,
                num_desc_components=num_desc_components,
                separator=separator,
                gbnum_in_front=gbnum_in_front,
                exclude_gbnum=exclude_gbnum)
    return taxon_set

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
        'pubmed',
        'protein',
        'nucleotide',
        'nuccore',
        'nucgss',
        'nucest',
        'structure',
        'genome',
        'biosystems',
        'books',
        'cancerchromosomes',
        'cdd',
        'gap',
        'dbvar',
        'domains',
        'epigenomics',
        'gene',
        'genomeprj',
        'gensat',
        'geo',
        'gds',
        'homologene',
        'journals',
        'mesh',
        'ncbisearch',
        'nlmcatalog',
        'omia',
        'omim',
        'pepdome',
        'pmc',
        'popset',
        'probe',
        'proteinclusters',
        'pcassay',
        'pccompound',
        'pcsubstance',
        'seqannot',
        'snp',
        'sra',
        'taxonomy',
        'toolkit',
        'toolkitall',
        'unigene',
        'unists',
        'linkoutpubmed',
        'linkoutseq',
        'linkoutother',
        ]

    class AccessionFetchError(Exception):

        def __init__(self, accession_ids):
            Exception.__init__(self, "Failed to retrieve accessions: %s" % (", ".join([str(s) for s in accession_ids])))

    def __init__(self,
            generate_labels=False,
            label_num_desc_components=3,
            label_separator='_',
            label_id_in_front=True,
            exclude_gbnum_from_label=False,
            sort_taxa_by_label=False):
        """
        *** DEPRECATED: use dendropy.interop.genbank.GenBankDna,
        dendropy.interop.genbank.GenBankRna, or
        dendropy.interop.genbank.GenBankProtein instead ***
        Instantiates a broker that queries NCBI and returns data.  If
        ``generate_labels`` is |True|, then appropriate labels for sequences
        will be automatically composed for each sequence based on the GenBank
        FASTA defline. ``label_num_desc_components`` specifies the number of
        components from the defline to use. ``label_separator`` specifies the
        string used to separate the different label components.
        ``label_id_in_front`` specifies whether the GenBank accession number
        should form the beginning (|True|) or tail (|False|) end of the
        label. ``sort_taxa_by_label`` specifies whether the sequences should be
        sorted by their final label values.
        """
        warnings.warn("This class (and parent module) has been deprecated. Use the classes provided in 'dendropy.interop.genbak' instead", DeprecationWarning)
        self.generate_labels = generate_labels
        self.label_num_desc_components = label_num_desc_components
        self.label_separator = label_separator
        self.label_id_in_front = label_id_in_front
        self.exclude_gbnum_from_label = exclude_gbnum_from_label
        self.sort_taxa_by_label = sort_taxa_by_label

    def fetch(self, db, ids, rettype):
        """
        Raw fetch. Returns file-like object opened for reading on string
        returned by query.
        """
        if textprocessing.is_str_type(ids):
            id_list = ids
        else:
            id_list = ",".join([str(i) for i in set(ids)])
        params = {'db': db,
                'id': id_list,
                'rettype': rettype,
                'retmode': 'text'}
        query_url = Entrez.BASE_URL + "/efetch.fcgi?" + urlio.urlencode(params)
        return urlio.read_url(query_url)

    def fetch_gbrecs_as_plaintext_dict(self,
            db,
            ids,
            key_field="accession"):
        db_name = "nucleotide"
        gb_recs_str = self.fetch(db=db, ids=ids, rettype="gb")
        gb_recs_str_list = re.split(r"^//$", gb_recs_str, flags=re.MULTILINE)
        gb_recs_str_list = [gb_rec for gb_rec in gb_recs_str_list
                    if gb_rec.replace("\n", "")]
        gb_recs_dict = {}
        for gb_str in gb_recs_str_list:
            if not gb_str:
                continue
            try:
                accession, gi = parse_accession_number_and_gi_from_gb(gb_str)
            except TypeError:
                print("---")
                print(gb_str)
                print("---")
                raise
            if key_field == "accession":
                gb_recs_dict[accession] = gb_str
            elif key_field == "gi":
                gb_recs_dict[gi] = gb_str
            else:
                raise NotImplementedError("Key field '%s' is not supported" % key_field)
        return gb_recs_dict

    def fetch_nucleotide_accessions(self,
            ids,
            prefix=None,
            verify=True,
            matrix_type=dendropy.DnaCharacterMatrix,
            **kwargs):
        """
        Returns a DnaCharacterMatrix object (or some other type, if specified
        by ``matrix_type`` argument) populated with sequences from the Entrez
        nucleotide database with accession numbers given by ``ids`` (a list of
        accession numbers). If ``prefix`` is given, it is pre-pended to all values
        given in the id list. Any other keyword arguments given are passed to
        the constructor of |DnaCharacterMatrix|.
        **Note that the order of records is *not* the same as the order of ids!!!**
        """
        if prefix is not None:
            id_list = ["%s%s" % (prefix,i) for i in ids]
        else:
            id_list = [str(i) for i in ids]
        results_str = self.fetch(db='nucleotide', ids=id_list, rettype='fasta')
        try:
            d = matrix_type.get_from_string(results_str, 'fasta', **kwargs)
        except DataParseError:
            sys.stderr.write("---\nNCBI Entrez Query returned:\n%s\n---\n" % results_str)
            raise
        for taxon in d.taxon_set:
            taxon.ncbi_defline = taxon.label
            taxon.ncbi_gi, taxon.ncbi_accession, taxon.ncbi_version = parse_ncbi_curation_info_from_defline(taxon.ncbi_defline)
        if verify:
            found_ids = set([t.ncbi_accession for t in d.taxon_set])
            missing_ids = set(id_list).difference(found_ids)
            found_ids = set([t.ncbi_gi for t in d.taxon_set])
            missing_ids = set(missing_ids).difference(found_ids)
            if len(missing_ids) > 0:
                raise Entrez.AccessionFetchError(missing_ids)
        if self.generate_labels:
            relabel_taxa_from_defline(d.taxon_set,
                    num_desc_components=self.label_num_desc_components,
                    separator=self.label_separator,
                    gbnum_in_front=self.label_id_in_front,
                    exclude_gbnum=self.exclude_gbnum_from_label)
        if self.sort_taxa_by_label:
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
        and, up to and *including* ``end``. If ``prefix`` is given, then it is
        pre-pended to the ids. Any other keyword arguments given are passed to
        thee constructor of |DnaCharacterMatrix|.
        """
        ids = range(first, last+1)
        return self.fetch_nucleotide_accessions(ids=ids, prefix=prefix, verify=verify, matrix_type=matrix_type, **kwargs)

