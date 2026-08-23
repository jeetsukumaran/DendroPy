#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
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

