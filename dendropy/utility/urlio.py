#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Communications using web/internet protocols.
"""

import sys
if sys.hexversion < 0x03000000:
    from urllib2 import Request
    from urllib2 import urlopen
    from urllib import urlencode
    from urllib2 import HTTPError
else:
    from urllib.request import Request
    from urllib.request import urlopen
    from urllib.parse import urlencode
    from urllib.error import HTTPError
import re
import locale

ENCODING = locale.getdefaultlocale()[1]

def read_url(url, strip_markup=False):
    """
    Return contents of url as string.
    """
    s = urlopen(url)
    text = s.read().decode(ENCODING)
    if strip_markup:
        return re.sub(r'<[^>]*?>', '', text)
    else:
        return text
