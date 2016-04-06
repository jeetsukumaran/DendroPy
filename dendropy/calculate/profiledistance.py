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

import collections
from dendropy.utility import constants
from dendropy.model import coalescent

class MeasurementProfile(object):

    def __init__(self, profile_data=None):
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

class TreeProfile(object):

    MEASUREMENT_NAMES = (
        "Edge.Lengths",
        "Patristic.Distances",
        "Patristic.Steps",
        "Node.Ages",
        "Coalescence.Intervals",
    )

    def __init__(self,
            tree,
            is_measure_node_ages=True,
            is_measure_coalescence_intervals=True,
            is_normalize=True,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            tree_phylogenetic_distance_matrix=None,
            tree_id=None,
            ):
        self.tree_id = tree_id
        self.measurement_profiles = collections.OrderedDict()

        if is_normalize:
            tree_length_nf = tree.length()
            if tree_length_nf == 0:
                tree_length_nf = 1.0
        else:
            tree_length_nf = 1.0
        self.measurement_profiles["Edge.Lengths"] = MeasurementProfile(
                profile_data=[float(e.length)/tree_length_nf for e in tree.postorder_edge_iter() if e.length is not None])
        if tree_phylogenetic_distance_matrix is None:
            tree_phylogenetic_distance_matrix = tree.phylogenetic_distance_matrix()
        self.measurement_profiles["Patristic.Distances"] = MeasurementProfile(
                profile_data=tree_phylogenetic_distance_matrix.distances(
                    is_weighted_edge_distances=True,
                    is_normalize_by_tree_size=is_normalize,))
        self.measurement_profiles["Patristic.Steps"] = MeasurementProfile(
                profile_data=tree_phylogenetic_distance_matrix.distances(
                    is_weighted_edge_distances=False,
                    is_normalize_by_tree_size=is_normalize,))
        if is_measure_node_ages:
            node_ages = tree.calc_node_ages(ultrametricity_precision=ultrametricity_precision)
            if is_normalize:
                s = sum(node_ages)
                node_ages = [a/s for a in node_ages]
            self.measurement_profiles["Node.Ages"] = MeasurementProfile(profile_data=node_ages,)
        if is_measure_coalescence_intervals:
            cf = coalescent.extract_coalescent_frames(
                    tree=tree,
                    ultrametricity_precision=ultrametricity_precision,)
            waiting_times = cf.values()
            if is_normalize:
                s = sum(waiting_times)
                waiting_times = [w/s for w in waiting_times]
            self.measurement_profiles["Coalescence.Intervals"] = MeasurementProfile(profile_data=node_ages,)

    def measure_distances(self, other_tree_profile):
        d = collections.OrderedDict()
        for pm_name in self.measurement_profiles:
            p1 = self.measurement_profiles[pm_name]
            p2 = other_tree_profile.measurement_profiles[pm_name]
            d[pm_name] = p1.euclidean_distance(p2)
        return d

class TreeProfileMatrix(object):

    def __init__(self,
            is_measure_node_ages=True,
            is_measure_coalescence_intervals=True,
            is_normalize=True,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            ):
        self.is_measure_node_ages = is_measure_node_ages
        self.is_measure_coalescence_intervals = is_measure_coalescence_intervals
        self.is_normalize = is_normalize
        self.ultrametricity_precision = ultrametricity_precision
        self.tree_profiles = []

    @property
    def measurement_names(self):
        names = [
            "Edge.Lengths",
            "Patristic.Distances",
            "Patristic.Steps",
        ]
        if self.is_measure_node_ages:
            names.append("Node.Ages")
        if self.is_measure_coalescence_intervals:
            names.append("Coalescence.Intervals")
        return names

    def add_tree(self,
            tree,
            tree_id=None,
            tree_phylogenetic_distance_matrix=None,
            ):
        if tree_id is None:
            tree_id = str(len(self.tree_profiles))
        profile = TreeProfile(
                tree=tree,
                tree_phylogenetic_distance_matrix=tree_phylogenetic_distance_matrix,
                tree_id=tree_id)
        self.tree_profiles.append(profile)
        return tree_id

    def compile(self):
        tree_profile_distances = collections.OrderedDict()
        for tree_profile_idx1, tree_profile1 in enumerate(self.tree_profiles[:-1]):
            tree_profile_distances[tree_profile1.tree_id] = collections.OrderedDict()
            for tree_profile_idx2, tree_profile2 in enumerate(self.tree_profiles[tree_profile_idx1+1:]):
                tree_profile_distances[tree_profile1.tree_id][tree_profile2.tree_id] = tree_profile1.measure_distances(tree_profile2)
        return tree_profile_distances
