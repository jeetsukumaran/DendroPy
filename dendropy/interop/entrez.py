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
Low-level wrappers around the NCBI E-Utilities. Primarily meant to open
file-like object handles on responses.
"""

import urllib

ENTREZ_BASE_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def entrez_get(**kwargs):
    query_url = ENTREZ_BASE_URL + "/efetch.fcgi?" + urllib.urlencode(kwargs)
    response = urllib.urlopen(query_url)
    return response

def entrez_post(**kwargs):
    pass

def efetch(db, ids, rettype, retmode="xml"):
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
    return entrez_get(**params)
