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

import math
import collections
from dendropy.utility import constants
from dendropy.model import coalescent

class MeasurementProfile(object):

    @staticmethod
    def _euclidean_distance(v1, v2, is_weight_values_by_comparison_profile_size=True):
        v1_size = len(v1)
        v2_size = len(v2)
        v1_idx = 0
        v2_idx = 0
        if v1_size > v2_size:
            v1_idx = v1_size - v2_size
            weight = float(v2_size)
        elif v2_size > v1_size:
            v2_idx = v2_size - v1_size
            weight = float(v1_size)
        else:
            weight = float(v1_size)
        if not is_weight_values_by_comparison_profile_size:
            weight = 1.0
        ss = 0.0
        while v1_idx < v1_size and v2_idx < v2_size:
            ss += pow(v1[v1_idx]/weight - v2[v2_idx]/weight, 2)
            v1_idx += 1
            v2_idx += 1
        return math.sqrt(ss)

    def __init__(self,
            profile_data=None,
            interpolation_method="piecewise_linear"):
        self.set_data(profile_data)
        self.fixed_size = 0
        self.interpolation_method = interpolation_method

    def add(self, value):
        self._profile_data.append(value)
        self._profile_data = sorted(self._profile_data)
        self._raw_data_size = len(self._profile_data)
        self._interpolated_profiles = {}

    def set_data(self, values):
        if values is None:
            values = []
        self._profile_data = sorted(values)
        self._raw_data_size = len(self._profile_data)
        self._interpolated_profiles = {}

    def __len__(self):
        return self._raw_data_size

    def distance(self,
            other,
            profile_size,
            is_weight_values_by_comparison_profile_size=True):
        if profile_size is None:
            profile_size = self._get_profile_comparison_size(other)
        v1 = self._get_profile_for_size(profile_size)
        v2 = other._get_profile_for_size(profile_size)
        return MeasurementProfile._euclidean_distance(v1, v2,
                is_weight_values_by_comparison_profile_size=is_weight_values_by_comparison_profile_size)
    # def get(self, idx):
    #     if idx >= self._raw_data_size:
    #         val = 0.0
    #     else:
    #         val = self._profile_data[idx]
    #         if val is None:
    #             val = 0.0
    #     return val

    def _get_profile_comparison_size(self, other):
        if self.fixed_size > 0 and other.fixed_size > 0:
            if self.fixed_size != other.fixed_size:
                raise ValueError("Comparing two profiles locked to different sizes: {} and {}".format(
                    self.fixed_size, other.fixed_size))
            return self.fixed_size
        elif self.fixed_size == 0 and other.fixed_size > 0:
            if other.fixed_size < self._raw_data_size:
                raise ValueError("Cannot interpolate points in current profile: current raw data size is {} but other profile is locked to a smaller fixed size ({})".format(
                    self._raw_data_size, other.fixed_size))
            return other.fixed_size
        elif self.fixed_size > 0 and other.fixed_size == 0:
            if self.fixed_size < other._raw_data_size:
                raise ValueError("Cannot interpolate points in other profile: other raw data size is {} but current profile is locked to a smaller fixed size ({})".format(
                    other._raw_data_size,
                    self.fixed_size,
                    ))
            return self.fixed_size
        else:
            return max(self._raw_data_size, other._raw_data_size)

    def _get_profile_for_size(self, profile_size):
        if not profile_size:
            return self._raw_data_size
        try:
            return self._interpolated_profiles[profile_size]
        except KeyError:
            self._interpolated_profiles[profile_size] = self._generate_interpolated_profile(profile_size)
            return self._interpolated_profiles[profile_size]

    def _generate_interpolated_profile(self, profile_size):
        if profile_size == self._raw_data_size:
            self._interpolated_profiles[profile_size] = list(self._profile_data)
            return self._interpolated_profiles[profile_size]
        if self._raw_data_size == 0:
            raise ValueError("No data in profile")
        if profile_size < self._raw_data_size:
            raise ValueError("Error interpolating points in profile: number of requested interpolated points ({}) is less than raw data size ({})".format(
                profile_size, self._raw_data_size))
        default_bin_size = int(profile_size / self._raw_data_size)
        if default_bin_size == 0:
            raise ValueError("Profile size ({}) is too small for raw data size ({}), resulting in a null bin size".
                    format(profile_size, self._raw_data_size))
        bin_sizes = [default_bin_size] * self._raw_data_size

        # // due to rounding error, default bin size may not be enough
        # // this hacky algorithm adds additional points to bins to make up the
        # // difference, trying to distribute the additional points along the
        # // length of the line
        diff = profile_size - (default_bin_size * self._raw_data_size)
        if diff > 0:
            dv = float(diff) / self._raw_data_size
            cv = 0.0
            for bin_idx in range(len(bin_sizes)):
                if diff <= 1.0:
                    break
                cv += dv
                if cv >= 1.0:
                    bin_sizes[bin_idx] += 1
                    diff -= 1.0
                    cv = cv - 1.0

        interpolated_profile = []
        if self.interpolation_method == "staircase":
            for original_data_value in self._profile_data:
                self._interpolate_flat(
                        interpolated_profile=interpolated_profile,
                        value=original_data_value,
                        num_points=default_bin_size)
        elif self.interpolation_method == "piecewise_linear":
            for bin_idx, original_data_value in enumerate(self._profile_data[:-1]):
                self._interpolate_linear(
                        interpolated_profile=interpolated_profile,
                        x1=bin_idx,
                        y1=self._profile_data[bin_idx],
                        y2=self._profile_data[bin_idx+1],
                        num_points=bin_sizes[bin_idx],
                        max_points=profile_size,
                        )
            interpolated_profile.append(self._profile_data[-1])
        return interpolated_profile

    def _interpolate_flat(self,
            interpolated_profile,
            value,
            num_points):
        interpolated_profile += [value] * num_points
        return interpolated_profile

    def _interpolate_linear(self,
            interpolated_profile,
            x1, y1, y2,
            num_points,
            max_points):
        assert num_points > 0
        slope = float(y2 - y1) / num_points
        for xi in range(num_points):
            if max_points and len(interpolated_profile) >= max_points:
                return interpolated_profile
            interpolated_profile.append((slope * xi) + y1)

class TreeProfile(object):

    def __init__(self,
            tree,
            is_measure_edge_lengths=True,
            is_measure_patristic_distances=False,
            is_measure_patristic_steps=False,
            is_measure_node_distances=True,
            is_measure_node_steps=True,
            is_measure_node_ages=True,
            is_measure_coalescence_intervals=True,
            is_normalize=True,
            ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            tree_phylogenetic_distance_matrix=None,
            tree_node_distance_matrix=None,
            tree_id=None,
            is_skip_normalization_on_zero_division_error=False,
            ):
        self.tree_id = tree_id
        self.is_measure_edge_lengths = is_measure_edge_lengths
        self.is_measure_patristic_distances = is_measure_patristic_distances
        self.is_measure_patristic_steps = is_measure_patristic_steps
        self.is_measure_node_distances = is_measure_node_distances
        self.is_measure_node_steps = is_measure_node_steps
        self.is_measure_node_ages = is_measure_node_ages
        self.is_measure_coalescence_intervals = is_measure_coalescence_intervals
        self.is_normalize = is_normalize
        self.is_skip_normalization_on_zero_division_error = is_skip_normalization_on_zero_division_error
        self.ultrametricity_precision = ultrametricity_precision
        self.measurement_profiles = collections.OrderedDict()
        self.compile(tree)

    def compile(self, tree,
            tree_phylogenetic_distance_matrix=None,
            tree_node_distance_matrix=None,
            ):
        if self.is_measure_edge_lengths:
            if self.is_normalize:
                tree_length_nf = tree.length()
                if tree_length_nf == 0:
                    tree_length_nf = 1.0
            else:
                tree_length_nf = 1.0
            self.measurement_profiles["Edge.Lengths"] = MeasurementProfile(
                    profile_data=[float(e.length)/tree_length_nf for e in tree.postorder_edge_iter() if e.length is not None])
        if self.is_measure_patristic_distances or self.is_measure_patristic_steps:
            if tree_phylogenetic_distance_matrix is None:
                tree_phylogenetic_distance_matrix = tree.phylogenetic_distance_matrix()
            if self.is_measure_patristic_distances:
                self.measurement_profiles["Patristic.Distances"] = MeasurementProfile(
                        profile_data=tree_phylogenetic_distance_matrix.distances(
                            is_weighted_edge_distances=True,
                            is_normalize_by_tree_size=self.is_normalize,))
            if self.is_measure_patristic_steps:
                self.measurement_profiles["Patristic.Steps"] = MeasurementProfile(
                        profile_data=tree_phylogenetic_distance_matrix.distances(
                        is_weighted_edge_distances=False,
                        is_normalize_by_tree_size=self.is_normalize,))
        if self.is_measure_node_distances or self.is_measure_node_steps:
            if tree_node_distance_matrix is None:
                tree_node_distance_matrix = tree.node_distance_matrix()
            if self.is_measure_node_distances:
                self.measurement_profiles["Node.Distances"] = MeasurementProfile(
                        profile_data=tree_node_distance_matrix.distances(
                            is_weighted_edge_distances=True,
                            is_normalize_by_tree_size=self.is_normalize,))
            if self.is_measure_node_steps:
                self.measurement_profiles["Node.Steps"] = MeasurementProfile(
                        profile_data=tree_node_distance_matrix.distances(
                        is_weighted_edge_distances=False,
                        is_normalize_by_tree_size=self.is_normalize,))
        if self.is_measure_node_ages:
            node_ages = tree.calc_node_ages(ultrametricity_precision=self.ultrametricity_precision)
            if self.is_normalize:
                s = sum(node_ages)
                try:
                    normalized_node_ages = [a/s for a in node_ages]
                    node_ages = normalized_node_ages
                except ZeroDivisionError as e:
                    if self.is_skip_normalization_on_zero_division_error:
                        pass
                    else:
                        raise
            self.measurement_profiles["Node.Ages"] = MeasurementProfile(profile_data=node_ages,)
        if self.is_measure_coalescence_intervals:
            cf = coalescent.extract_coalescent_frames(
                    tree=tree,
                    ultrametricity_precision=self.ultrametricity_precision,)
            waiting_times = cf.values()
            if self.is_normalize:
                s = sum(waiting_times)
                try:
                    normalized_waiting_times = [w/s for w in waiting_times]
                    waiting_times = normalized_waiting_times
                except ZeroDivisionError as e:
                    if self.is_skip_normalization_on_zero_division_error:
                        pass
                    else:
                        raise
            self.measurement_profiles["Coalescence.Intervals"] = MeasurementProfile(profile_data=waiting_times,)

    @property
    def measurement_names(self):
        names = []
        if self.is_measure_edge_lengths:
            names.append("Edge.Lengths")
        if self.is_measure_patristic_distances:
            names.append("Patristic.Distances")
        if self.is_measure_patristic_steps:
            names.append("Patristic.Steps")
        if self.is_measure_node_distances:
            names.append("Node.Distances")
        if self.is_measure_node_steps:
            names.append("Node.Steps")
        if self.is_measure_node_ages:
            names.append("Node.Ages")
        if self.is_measure_coalescence_intervals:
            names.append("Coalescence.Intervals")
        return names

    def measure_distances(self, other_tree_profile,
            profile_size=None,
            is_weight_values_by_comparison_profile_size=True):
        d = collections.OrderedDict()
        for pm_name in self.measurement_profiles:
            p1 = self.measurement_profiles[pm_name]
            p2 = other_tree_profile.measurement_profiles[pm_name]
            d[pm_name] = p1.distance(p2,
                    profile_size=None,
                    is_weight_values_by_comparison_profile_size=is_weight_values_by_comparison_profile_size)
        return d

