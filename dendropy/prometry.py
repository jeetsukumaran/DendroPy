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
Profile distances.
"""

class Profile(object):

    def __init__(self,
            key=None,
            profile_data=None):
        self.key = key
        if profile_data is None:
            self._profile_data = []
        else:
            self._profile_data = sorted(profile_data)

    def add(self, value):
        self._profile_data.append(value)
        self._profile_data = sorted(self._profile_data)

    def __len__(self):
        return len(self._profile_data)

    def get(self, idx):
        if idx >= len(self._profile_data):
            val = 0.0
        else:
            val = self._profile_data[idx]
            if val is None:
                val = 0.0
        return val

    def euclidean_distance(self, other):
        distance = 0.0
        for i in range(max(len(self._profile_data), len(other._profile_data))):
            di = pow(self.get(i) - other.get(i), 2)
            distance += di
        return distance

class ProfileMatrix(object):

    def __init__(self, title):
        self.title = title
        self._profiles = []

    def add(self, key, profile_data):
        profile = Profile(key=key,
                profile_data=profile_data)
        self._profiles.append(profile)

    def __len__(self):
        return len(self._profiles)

    def _calc_distances(self, dist_func):
        results = {}
        for idx1, profile1 in enumerate(self._profiles[:-1]):
            for idx2, profile2 in enumerate(self._profiles[idx1+1:]):
                d = dist_func(profile1, profile2)
                results[frozenset([profile1, profile2])] = d
        return results

    def calc(self):
        results = Distances(title=self.title)
        euclidean_dist_func = lambda x, y: x.euclidean_distance(y)
        d = self._calc_distances(euclidean_dist_func)
        results.add("Euclidean.Raw", d)
        return results

class Distances(object):

    def __init__(self, title):
        self.title = title
        self.distance_names = []
        self.distances = {}
        self.comparisons = set()

    def add(self, name, distance_dict):
        """
            - `name`
                string describing distance type (e.g., "Raw.Euclidean")
            - `distance_dict`
                dictionary with frozen set of Profile objects being compared
                as keys, and distance value as values
        """
        self.distance_names.append(name)
        self.distances[name] = distance_dict
        self.comparisons.update(distance_dict.keys())

    def get(self, distance_name, comparison):
        return self.distances[distance_name][frozenset(comparison)]


def summarize_profile_matrices(profile_matrices,
        out,
        include_header_row=True):
    distance_groups = {}
    distance_group_titles = []
    distance_names = set()
    comparisons_set = set()
    for pm in profile_matrices:
        d = pm.calc()
        distance_groups[d.title] = d
        distance_group_titles.append(d.title)
        distance_names.update(d.distance_names)
        comparisons_set.update(d.comparisons)
    distance_group_titles = sorted(distance_group_titles)
    distance_names = sorted(distance_names)
    comparisons_list = []
    for comparison in comparisons_set:
        comparisons_list.append( sorted(comparison, key=lambda x: x.key) )
    comparisons_list = sorted(comparisons_list, key=lambda x: x[0].key + x[1].key)
    if include_header_row:
        row_parts = ["P1", "P2"]
        for dist_title in distance_group_titles:
            for dist_name in distance_names:
                row_parts.append("{}.{}".format(dist_title, dist_name))
        out.write("\t".join(row_parts))
        out.write("\n")
    for comparisons in comparisons_list:
        row_parts = [comparisons[0].key, comparisons[1].key]
        for dist_title in distance_group_titles:
            for dist_name in distance_names:
                distance = distance_groups[dist_title].get(dist_name, comparisons)
                row_parts.append("{}".format(distance))
        out.write("\t".join(row_parts))
        out.write("\n")
