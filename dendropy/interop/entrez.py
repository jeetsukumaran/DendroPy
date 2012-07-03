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

ENTREZ_EUTILS_BASE_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def efetch(db, ids, rettype, retmode="xml", email=None):
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
    if email is not None:
        params["email"] = email
    query_url = ENTREZ_EUTILS_BASE_URL + "/efetch.fcgi?" + urllib.urlencode(params)
    response = urllib.urlopen(query_url)
    return response

def get_taxonomy(**kwargs):
    params = dict(kwargs)
    params["db"] = "taxonomy"
    # if "rettype" not in params:
    #     params["rettype"] = "gbc"
    if "retmode" not in params:
        params["retmode"] = "xml"
    if "email" in kwargs and kwargs["email"] is None:
        del params["email"]
    query_url = "http://www.ncbi.nlm.nih.gov/sites/entrez?" + urllib.urlencode(params)
    response = urllib.urlopen(query_url)
    return response

# # def efetch(db, ids, rettype, retmode="xml", email=None):
# def efetch(db, **kwargs):
#     """
#     Raw fetch. Returns file-like object opened for reading on string
#     returned by query.
#     """
#     params = dict(kwargs)
#     params["db"] = db
#     if "email" in kwargs and kwargs["email"] is None:
#         del params["email"]
#     if "id" in params:
#         params["ids"] = params["id"]
#         del params["id"]
#     if "rettype" not in params:
#         params["rettype"] = "gbc"
#     if "retmode" not in params:
#         params["retmode"] = "xml"
#     if "ids" in params:
#         ids = params["ids"]
#         if isinstance(ids, str):
#             id_list = ids
#         else:
#             id_list = ",".join([str(i) for i in list(ids)])
#         del params["ids"]
#         params["id"] = id_list
#     return entrez_get(**params)
