#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Communications using web/internet protocols.
"""

from dendropy.utility import textprocessing
from dendropy.utility import error

from urllib.request import urlopen
import re

def read_url(url, strip_markup=False):
    """
    Return contents of url as string.
    """
    s = urlopen(url)
    text = textprocessing.bytes_to_text(s.read())
    if strip_markup:
        return re.sub(r'<[^>]*?>', '', text)
    else:
        return text

def post_request(url, data=None, **kwargs):
    try:
        import requests
    except ImportError:
        msg = ("\n"
              "This operation requires installation of the 'Requests' library:\n"
              "\n"
              "    http://docs.python-requests.org/en/master/ \n"
              "\n")
        raise error.LibraryDependencyError(msg)
    response = requests.post(url=url, data=data, **kwargs)
    return response

