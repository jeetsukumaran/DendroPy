#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Wrappers for interacting with the APE library for R.
"""

import tempfile
import re
import sys
import dendropy

try:
    from rpy2 import robjects
    from rpy2.rinterface import RRuntimeError
    _R = robjects.r
    _R('library(ape)')
except ImportError:
    sys.stderr.write("rpy2 not installed: APE interoperability not available")
except RRuntimeError:
    sys.stderr.write("APE library not installed: APE interoperability not available")
else:

    def as_ape_object(o):
        """
        Returns `o` as an ape object.
        """
        kwargs = {}
        if isinstance(o, dendropy.TreeList):
            kwargs['keep.multi'] = True
            text = o.as_string("newick")
            return _R['read.tree'](text=text, **kwargs)
        elif isinstance(o, dendropy.Tree):
            kwargs['keep.multi'] = False
            text = o.as_string("newick")
            return _R['read.tree'](text=text, **kwargs)
        elif isinstance(o, dendropy.CharacterMatrix):
            f = tempfile.NamedTemporaryFile()
            o.write_to_stream(f, "nexus", simple=True)
            f.flush()
            return _R['read.nexus.data'](f.name)
        else:
            return robjects.default_py2ri(o)

    def as_r_vector(o, val_type):
        if isinstance(o, dict):
            keys = o.keys()
            vals = [o[k] for k in keys]
            robj = as_r_vector(vals, val_type=val_type)
            robj.setnames(keys)
        else:
            if val_type == int:
                robj = robjects.IntVector(o)
            elif val_type == float:
                robj = robjects.FloatVector(o)
            elif val_type == bool:
                robj = robjects.BoolVector(o)
            else:
                robj = robjects.RVector(o)
        return robj

    def as_dendropy_object(o, taxon_set=None):
        """
        Returns a DendroPy object corresponding to the ape object `o`. If `o` is
        a single tree (i.e., `phylo`), then a DendroPy Tree is returned. If `o` is
        a list of trees (i.e., a `multiPhylo` object, or list of `phylo` objects),
        then a DendroPy TreeList is returned.
        """
        if o.rclass[0] == "multiPhylo":
            f = tempfile.NamedTemporaryFile()
            _R['write.nexus'](o, file=f.name)
            return dendropy.TreeList.get_from_path(f.name, "nexus", taxon_set=taxon_set)
        elif o.rclass[0] == "phylo":
            f = tempfile.NamedTemporaryFile()
            _R['write.nexus'](o, file=f.name)
            return dendropy.Tree.get_from_path(f.name, "nexus", taxon_set=taxon_set)
        elif o.rclass[0] == "list":
            f = tempfile.NamedTemporaryFile()
            _R['write.nexus.data'](o, file=f.name)
    #        print open(f.name, "r").read()
            d = dendropy.DataSet.get_from_path(f.name, "nexus", taxon_set=taxon_set)
            if len(d.char_matrices) == 0:
                raise ValueError("No character data found")
            elif len(d.char_matrices) == 1:
                return d.char_matrices[0]
            else:
                raise ValueError("Multiple character matrices returned")
        else:
            return robjects.default_ri2py(o)

    def exec_and_capture(rfunc, *args, **kwargs):
        stdoutf = tempfile.NamedTemporaryFile()
        stderrf = tempfile.NamedTemporaryFile()
        _R('sink("%s")' % stdoutf.name)
        _R('zz = file("%s", "wt")' % stderrf.name)
        _R('sink(zz, type="message")')
        rfunc(*args, **kwargs)
        _R('sink(type="message")')
        _R('sink()')
        i = open(stdoutf.name, "rU")
        stdout = i.read()
        i = open(stderrf.name, "rU")
        stderr = i.read()
        return stdout, stderr

    def bd_ext(t, num_species_node_attr='num_species'):
        """
        This function fits by maximum likelihood a birth-death model to
        the combined phylogenetic and taxonomic data of a given clade. The
        phylogenetic data are given by a tree, `t`, and the taxonomic data by
        an attribute `num_species` of each of the taxa in the tree.
        Returns dictionary, where keys are param names and values are param
        values.
        """
        taxon_num_species = []
        for taxon in t.taxon_set:
            taxon_num_species.append(taxon.num_species)
    #    taxon_num_species_map = {}
    #    for taxon in t.taxon_set:
    #        taxon_num_species_map[taxon.label.replace(" ", "_")] = taxon.num_species
    #    taxon_num_species_map = [10, 47, 69, 214, 161, 17,
    #            355, 51, 56, 10, 39, 152,
    #            6, 143, 358, 103, 319,
    #            23, 291, 313, 196, 1027, 5712]
        stdout, stderr = exec_and_capture(_R['bd.ext'], as_ape_object(t), as_r_vector(taxon_num_species, int))
        patterns = {
            'deviance' : '\s*Deviance: ([\d\-\.Ee\+]+).*',
            'log-likelihood' : '\s*Log-likelihood: ([\d\-\.Ee\+]+)',
            'd/b' : '\s*d / b = ([\d\-\.Ee\+]+)',
            'd/b s.e.' : '\s*d / b = .* StdErr = ([\d\-\.Ee\+]+)',
            'b-d' : '\s*b - d = ([\d\-\.Ee\+]+)',
            'b-d s.e.' : '\s*b - d = .* StdErr = ([\d\-\.Ee\+]+)',
        }
        results = {}
        for k, v in patterns.items():
            m = re.findall(v, stdout)
            if m:
                results[k] = float(m[0])
        return results

    def birthdeath(t):
        """
        This function fits by maximum likelihood a birth-death model to
        the branching times computed from a phylogenetic tree using the
        method of Nee et al. (1994).
        Returns dictionary, where keys are param names and values are param
        values.
        """
        _R('options(warn=-99)')
        rval = _R['birthdeath'](as_ape_object(t))
        _R('options(warn=0)')
        results = {}
        names = [n for n in rval.names]
        results['deviance'] = float(rval[names.index('dev')][0])
        results['log-likelihood'] = results['deviance'] / -2
        para = rval[names.index('para')]
        para_names = [n for n in para.names]
        results['d/b'] = float(para[para_names.index('d/b')])
        results['b-d'] = float(para[para_names.index('b-d')])
        se = rval[names.index('se')]
        se_names = [n for n in se.names]
        results['d/b s.e.'] = float(se[se_names.index('d/b')])
        results['b-d s.e.'] = float(se[se_names.index('b-d')])
        return results


