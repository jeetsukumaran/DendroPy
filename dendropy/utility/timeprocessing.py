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

import datetime

def pretty_timestamp(t=None, style=0):
    if t is None:
        t = datetime.datetime.now()
    if style == 0:
        return t.strftime("%Y-%m-%d")
    elif style == 1:
        return t.strftime("%Y%m%d%H%M%S")
    else:
        return t.strftime("%Y%m%d%H%M%S%f")

def pretty_elapsed_datetime(t, fill=False):
    parts = []
    _render = lambda f, value: "{} {}{}".format(value, f, "" if value == 1 else "s")
    if t.day or fill:
        parts.append(_render("day", t.day))
    if t.hour or fill:
        parts.append(_render("hour", t.hour))
    if t.minute or fill:
        parts.append(_render("minute", t.minute))
    secs = t.second + float(t.microsecond)/1000000
    if secs or fill:
        s = _render("second", secs)
        if parts:
            parts.append("and {}".format(s))
        else:
            parts.append(s)
    return ", ".join(parts)

def parse_timedelta(td):
    hours = (td.days * 24) + td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    # seconds = ((td.seconds % 3600) % 60) + float(td.microseconds)/1000000
    seconds = ((td.seconds % 3600) % 60) + float(td.microseconds)/1000000
    return hours, minutes, seconds

def pretty_timedelta(td, fill=False):
    hours, minutes, seconds = parse_timedelta(td)
    parts = []
    _render = lambda f, value: "{} {}{}".format(value, f, "" if value == 1 else "s")
    if hours or fill:
        parts.append(_render("hour", hours))
    if minutes or fill:
        parts.append(_render("minute", minutes))
    if seconds or fill:
        s = _render("second", seconds)
        if parts:
            parts.append("and {}".format(s))
        else:
            parts.append(s)
    if not parts:
        return "0 seconds"
    return ", ".join(parts)
