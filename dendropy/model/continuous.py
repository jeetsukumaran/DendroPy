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

"""
Models and modeling of continuous character evolution.
"""

import math
from functools import reduce
import operator
import dendropy
from dendropy.utility import GLOBAL_RNG

class PhylogeneticIndependentConstrasts(object):
    """
    Phylogenetic Independent Contrasts.

    References:

        -   Felsenstein, J. 1985. Phylogenies and the comparative method. American
            Naturalist 125:1-15.
        -   Garland, T., Jr., Jr., A. F. Bennett, and E. L. Rezende. 2005.
            Phylogenetic approaches in comparative physiology. Journal of
            Experimental Biology 208:3015-3035.

    """

    def __init__(self,
            tree,
            char_matrix,
            polytomy_strategy=None):
        """

        Parameters
        ----------
        tree : |Tree| object
            Tree to use.
        char_matrix : |ContinuousCharacterMatrix|
            ContinuousCharacterMatrix that is the source of the data
        polytomy_strategy
            One of: 'error', 'ignore', 'resolve'.

                'error'
                    Throws an error if tree has polytomies.
                'ignore'
                    No error, but raw contrasts will not be calculated for
                    polytomies.
                'resolve'
                    Randomly resolve polytomies.

            Defaults to 'error' if not specified or set to None.
        """
        self._tree = None
        self._char_matrix = None
        self._is_dirty = None
        self._is_fully_analyzed = False
        self._polytomy_strategy = None
        self._character_contrasts = {}
        self._set_polytomy_strategy(polytomy_strategy)
        self.tree = tree
        self.char_matrix = char_matrix

    def _get_polytomy_strategy(self):
        return self._polytomy_strategy
    def _set_polytomy_strategy(self, polytomy_strategy):
        if polytomy_strategy is None:
            self._polytomy_strategy = 'error'
        else:
            polytomy_strategy = polytomy_strategy.lower()
            if polytomy_strategy not in ['error', 'ignore', 'resolve']:
                raise ValueError("Invalid polytomy strategy: '%s'" % polytomy_strategy)
            else:
                self._polytomy_strategy = polytomy_strategy
    polytomy_strategy = property(_get_polytomy_strategy, None)

    def _get_tree(self):
        if not self._is_fully_analyzed:
            analyzed_chars = self._character_contrasts.keys()
            for idx in range(len(self.char_matrix[0])):
                if idx in analyzed_chars:
                    continue
                self._get_contrasts(idx)
            self._is_fully_analyzed = True
        return self._tree
    def _set_tree(self, tree):
        self._tree = dendropy.Tree(tree)
        if self._polytomy_strategy == 'resolve':
            self._tree.resolve_polytomies()
        self.is_dirty = True
    tree = property(_get_tree, _set_tree, None, """\
            This tree will have an attribute added to each node, ``pic``. This
            attribute will be a dictionary with character (column) index as
            keys. Each column index will map to another dictionary that has the
            following keys (and values):

                - ``pic_state_value``
                - ``pic_state_variance``
                - ``pic_contrast_raw``
                - ``pic_contrast_variance``
                - ``pic_contrast_standardized``
                - ``pic_edge_length_error``
                - ``pic_corrected_edge_length``
            """)

    def _get_char_matrix(self):
        return self._char_matrix
    def _set_char_matrix(self, char_matrix):
        self._char_matrix = char_matrix
        self.is_dirty = True
    char_matrix = property(_get_char_matrix, _set_char_matrix)

    def _get_is_dirty(self):
        return self._is_dirty
    def _set_is_dirty(self, is_dirty):
        self._is_dirty = is_dirty
        if self._is_dirty:
            self._character_contrasts = {}
            self._is_fully_analyzed = False
    is_dirty = property(_get_is_dirty, _set_is_dirty)

    def _get_contrasts(self, character_index):
        """
        Main work-horse method. If needed, adds an entry to
        self._character_constrants, with key being the character index, and a
        value being another dictionary that contains the constrast information.
        This second dictionary has the node's id as a key and as a value the a
        dictionary with the following:

                - ``pic_state_value``
                - ``pic_state_variance``
                - ``pic_contrast_raw``
                - ``pic_contrast_variance``
                - ``pic_contrast_standardized``
                - ``pic_edge_length_error``
                - ``pic_corrected_edge_length``

        """
        if character_index in self._character_contrasts:
            return self._character_contrasts[character_index]
        all_results = {}
        for nd in self._tree.postorder_node_iter():
            nd_results = {}
            child_nodes = nd.child_nodes()
            if len(child_nodes) == 0:
                nd_results['pic_state_value'] = self._char_matrix[nd.taxon][character_index]
                nd_results['pic_state_variance'] = None
                nd_results['pic_contrast_raw'] = None
                nd_results['pic_contrast_variance'] = None
                nd_results['pic_contrast_standardized'] = None
                nd_results['pic_edge_length_error'] = 0.0
                nd_results['pic_corrected_edge_length'] = nd.edge.length
            elif len(child_nodes) == 1:
                # root node?
                nd_results['pic_state_value'] = None
                nd_results['pic_state_variance'] = None
                nd_results['pic_contrast_raw'] = None
                nd_results['pic_contrast_variance'] = None
                nd_results['pic_contrast_standardized'] = None
                nd_results['pic_edge_length_error'] = None
                nd_results['pic_corrected_edge_length'] = None
            else:
                state_vals = []
                corrected_edge_lens = []
                actual_edge_lens = []
                for cnd in child_nodes:
                    state_vals.append(all_results[cnd._track_id]['pic_state_value'])
                    actual_edge_lens.append(cnd.edge.length)
                    if all_results[cnd._track_id]['pic_corrected_edge_length'] is not None:
                        corrected_edge_lens.append(all_results[cnd._track_id]['pic_corrected_edge_length'])
                    else:
                        corrected_edge_lens.append(cnd.edge.length)
                n = len(state_vals)
                numerator_fn = lambda i : (1.0/corrected_edge_lens[i]) * state_vals[i]
                denominator_fn = lambda i  : 1.0/corrected_edge_lens[i]
                nd_results['pic_state_value'] = \
                        sum(numerator_fn(i) for i in range(n)) \
                        / sum(denominator_fn(i) for i in range(n))

                sum_of_child_edges = sum(corrected_edge_lens)
                prod_of_child_edges = reduce(operator.mul, corrected_edge_lens)
                nd_results['pic_edge_length_error'] = (  prod_of_child_edges / (sum_of_child_edges) )
                if nd.edge.length is not None:
                    nd_results['pic_corrected_edge_length'] = nd.edge.length + nd_results['pic_edge_length_error']
                else:
                    nd_results['pic_corrected_edge_length'] = None
                nd_results['pic_state_variance'] = nd_results['pic_corrected_edge_length']

                if len(child_nodes) != 2:
                    if self._polytomy_strategy == "ignore":
                        nd_results['pic_contrast_raw'] = None
                        nd_results['pic_contrast_standardized'] = None
                        nd_results['pic_contrast_variance'] = sum_of_child_edges
                    else:
                        raise ValueError("Tree is not fully-bifurcating")
                else:
                    nd_results['pic_contrast_raw'] = state_vals[0] - state_vals[1]
                    nd_results['pic_contrast_standardized'] = nd_results['pic_contrast_raw'] / (sum_of_child_edges ** 0.5)
                    nd_results['pic_contrast_variance'] = sum_of_child_edges

            nd._track_id = id(nd) # will get cloned
            all_results[nd._track_id] = nd_results
            try:
                nd.pic[character_index] = dict(nd_results)
            except AttributeError:
                nd.pic = {character_index: dict(nd_results)}
        self._character_contrasts[character_index] = dict(all_results)
        return self._character_contrasts[character_index]

    def contrasts_tree(self,
            character_index,
            annotate_pic_statistics=True,
            state_values_as_node_labels=False,
            corrected_edge_lengths=False):
        """
        Returns a Tree object annotated with the following attributes added
        to each node (as annotations to be serialized if
        ``annotate_pic_statistics`` is True):

            - ``pic_state_value``
            - ``pic_state_variance``
            - ``pic_contrast_raw``
            - ``pic_contrast_variance``
            - ``pic_contrast_standardized``
            - ``pic_edge_length_error``
            - ``pic_corrected_edge_length``

        """
        contrasts = self._get_contrasts(character_index)
        tree = dendropy.Tree(self._tree)
        for nd in tree.postorder_node_iter():
            nd_results = contrasts[nd._track_id]
            for k, v in nd_results.items():
                setattr(nd, k, v)
                if annotate_pic_statistics:
                    nd.annotations.add_bound_attribute(k)
            if corrected_edge_lengths and nd_results['pic_corrected_edge_length'] is not None:
                nd.edge.length = nd_results['pic_corrected_edge_length']
            if state_values_as_node_labels:
                nd.label = str(nd_results['pic_state_value'])
        return tree

def evolve_continuous_char(node, rng=None, **kwargs):
    """
    Takes a node and a random number generator object, ``rng`` This function
    "evolves" a set of rates on the subtree descending from the  ``node``.

    kwargs keys that are used are:

    ``roeotroe``
        the rate of evolution of the rate of evolution. This
        parameter that controls the degree of deviations away from the
        molecular clock.
    ``min_rate``
        is the minimum rate (default None)
    ``max_rate``
        is the maximum rate (default None),
    ``model``
        is a string specifying the name of the model. Currently only the
        KTB (Kishino, Thorne, Bruno) is supported
    ``time_attr``
        is a string that specifies the name of the attribute
        that returns the branch length in terms of time for a node. The
        default is "edge_length"
    ``val_attr``
        is the string that specifies the name of the attribute
        used to hold the value that is evolving along the nodes.  The root
        of the subtree is assumed to have this field on calling of the
        function.  On success all nodes in the subtree will have the
        attribute.  The default is "mutation_rate"
    ``mean_val_attr``
        if specified this is string that gives the name of
        attribute in each node that is mean value for the branch (default is
        None). This is filled in after time_attr and val_attr are read,
        so it is permissible to have this attribute match one of thos
        strings (although it will make the model odd if the mean_val_attr
        is the same as the val_attr)
        ``constrain_rate_mode``
        controls the behavior when the minimum or maximum rate is
        simulated. The choices are "crop", and "linear_bounce" "crop" means
        that the rate is set to the most extreme value allowed.
        "linear_bounce" considers the path of evolution of rate to be a
        simple line from the ancestor's rate to the proposed rate. The
        point at which the path crosses the extreme value is determined and
        the rate is "reflected off" the limiting rate at that point.  This
        causes the rate to avoid the extreme values more than a simulation
        of small time slices that simply rejects illegal rates.

    Currently the only model supported is the one of Kishino, Thorne, and Bruno.
    "Performance of a Divergence Time Estimation Method under a Probabilistic
    Model of Rate Evolution." Molecular Biology and Evolution (2001) vol. 18
    (3) pp. 352-361. This model is specified by the code "KTB". A node's rate
    is a log-normal variate with variance determined by the product of the
    duration of the branch and the roeotroe parameter.  The mean of the
    distribution is chosen such that mean of the log-normal is identical to the
    rate at the parent. The mean_rate for the branch is the average of the rates
    at the endpoints.

    """
    if rng is None:
        rng = GLOBAL_RNG
    nd_iter = node.preorder_iter()
    # skip the first node -- it should already have a rate
    next(nd_iter)
    if kwargs.get("model", "KTB").upper() != "KTB":
        raise ValueError("Only the Kishino-Thorne-Bruno model is supported at this time")
    val_attr = kwargs.get("val_attr", "mutation_rate")
    if not val_attr:
        raise ValueError("val_attr cannot be an empty string")
    time_attr = kwargs.get("time_attr", "edge_length")
    mean_val_attr = kwargs.get("mean_val_attr")
    constrain_rate_mode = kwargs.get("constrain_rate_mode", "crop").lower()
    if constrain_rate_mode not in ["crop", "linear_bounce"]:
        raise ValueError('Only "crop" and "linear_bounce" are supported at this time')
    roeotroe = kwargs.get("roeotroe", 1.0)
    min_rate = kwargs.get("min_rate", 0.0)
    if min_rate < 0.0:
        raise ValueError("min_rate cannot be less than 0")
    max_rate = kwargs.get("max_rate")
    anc_rate = getattr(node, val_attr)
    if max_rate is not None:
        if min_rate is not None:
            if min_rate > max_rate:
                raise ValueError("max_rate must be greater than the min_rate")
            if min_rate == max_rate:
                for nd in nd_iter:
                    setattr(nd, val_attr, min_rate)
                    if mean_val_attr:
                        # here we assume that the rate changed from the
                        #   ancestral rate to the only allowed rate
                        #   instantaneously, so the mean rat is min_rate
                        setattr(nd, mean_val_attr, min_rate)
                return
        if max_rate <= 0.0:
            raise ValueError("max_rate must be positive")
        if anc_rate > max_rate:
            raise ValueError("rate for the incoming node is > max_rate")
    if (min_rate is not None) and anc_rate < min_rate:
        raise ValueError("rate for the incoming node is > max_rate")

    if constrain_rate_mode == "crop":
        rate_fn = _calc_KTB_rates_crop
    else:
        rate_fn = _calc_KTB_rates_linear_bounce
    for nd in nd_iter:
        starting_rate = getattr(nd.parent_node, val_attr)
        duration = getattr(nd, time_attr)
        r, mr  = rate_fn(starting_rate, duration, roeotroe, rng, min_rate, max_rate)
        setattr(nd, val_attr, r)
        if mean_val_attr:
            setattr(nd, mean_val_attr, mr)

def _bounce_constrain(start_x, x, min_x=None, max_x=None):
    """Returns the value of variable and its mean value over a path.
    We assume that some variable started at ``start_x`` and moved toward ``x``, but
    has to bounce of barriers specified by ``min_x`` and ``max_x``.

    ``x`` determines the direction and magnitude of the change.

    ``start_x`` must fall in the legal range (between the min and max). If
    ``x`` is also legal, then (x, (x + start_x)/2.0) will be returned reflecting
    the fact that the arithmetic mean of the endpoints represents the mean value
    of the variable if it took a direct path (at constant rate).
    """

    if max_x is not None and min_x is not None:
        assert(max_x > min_x)
    gt_max = (max_x is not None  and x > max_x)
    lt_min = (min_x is not None and x < min_x)
    prev_x = start_x
    prop_dur_remaining = 1.0
    mx = 0.0
    while gt_max or lt_min:
        if gt_max:
            p_changing = (max_x - prev_x)/(x - prev_x)
            mean_changing = (prev_x + max_x)/2.0
            mx += p_changing*prop_dur_remaining*mean_changing
            prop_dur_remaining *= 1.0 - p_changing
            x = 2*max_x - x
            lt_min = (min_x is not None and x < min_x)
            prev_x =  max_x
            gt_max = False
        if lt_min:
            p_changing = (prev_x - min_x)/(prev_x - x)
            mean_changing = (prev_x + min_x)/2.0
            mx += prop_dur_remaining*p_changing*mean_changing
            prop_dur_remaining *= 1.0 - p_changing
            x = 2*min_x - x
            lt_min = False
            gt_max = (max_x is not None  and x > max_x)
            prev_x =  min_x
    mean_changing = (prev_x + x)/2.0
    mx += mean_changing*prop_dur_remaining
    return x, mx

def _calc_TKP_rate(starting_rate, duration, roeotroe, rng):
    """
    Returns a simulated rate for the head node of a tree when:

        * the tail node has rate ``starting_rate``
        * the time duration of the edge is ``duration``
        * the rate of evolution of the rate of evolution is ``roeotroe`` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    ``rng`` is a random number generator.

    The model used to generate the rate is the one described by Thorne, Kishino,
    and Painter 1998.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of ``starting_rate``
    The variance of the normal distribution for the logarithm of the ending rate
        is the product of ``duration`` and ``roeotroe``
    """
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        mu = math.log(starting_rate)
        return rng.lognormvariate(mu, math.sqrt(rate_var))
    return starting_rate

def _calc_KTB_rate(starting_rate, duration, roeotroe, rng):
    """
    Returns a simulated rate for the head node of a tree when:
        * the tail node has rate ``starting_rate``
        * the time duration of the edge is ``duration``
        * the rate of evolution of the rate of evolution is ``roeotroe`` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    ``rng`` is a random number generator.

    The model used to generate the rate is the one described by Kishino, Thorne,
    and Bruno 2001.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of ``starting_rate``
    The variance of the normal distribution for the logarithm of the ending rate
        is the product of ``duration`` and ``roeotroe``
    """
    if starting_rate <= 0.0:
        raise ValueError("starting_rate must be positive in the KTB model")
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        # Kishino, Thorne and Bruno corrected the tendency for the rate to
        #   increase seen in teh TKP, 1998 model
        mu = math.log(starting_rate) - (rate_var/2.0)
        return rng.lognormvariate(mu, math.sqrt(rate_var))
    return starting_rate

def _calc_KTB_rates_crop(starting_rate, duration, roeotroe, rng,  min_rate=None, max_rate=None):
    """Returns a descendant rate and mean rate according to the Kishino, Thorne,
    Bruno model.  Assumes that the min_rate <= starting_rate <= max_rate if a max
    and min are provided.
    rate is kept within in the [min_rate, max_rate] range by cropping at these
    values and acting is if the cropping occurred at an appropriate point
    in the duration of the branch (based on a linear change in rate from the
    beginning of the random_variate drawn for the end of the branch).
    """
    if roeotroe*duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError("Parent rate is out of bounds, but no rate change is possible")
    r = _calc_KTB_rate(starting_rate, duration, roeotroe, rng)
    if max_rate and r > max_rate:
        assert(starting_rate <= max_rate)
        p_changing =  (max_rate - starting_rate)/(r - starting_rate)
        mean_changing = (starting_rate + max_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*max_rate
        return max_rate, mr
    elif min_rate and r < min_rate:
        assert(starting_rate >= min_rate)
        p_changing = (starting_rate - min_rate)/(starting_rate - r)
        mean_changing = (starting_rate + min_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*min_rate
        return min_rate, mr
    return r, (starting_rate + r)/2.0

def _calc_KTB_rates_linear_bounce(starting_rate, duration, roeotroe, rng,  min_rate=0.0, max_rate=None):
    """Returns a descendant rate and mean rate according to the Kishino, Thorne,
    Bruno model.  Assumes that the min_rate <= starting_rate <= max_rate if a max
    and min are provided.
    The rate is kept within in the [min_rate, max_rate] range by "bouncing" off
    of the barriers, where the "collision" is estimated by assuming a linear
    change in rate from the beginning of the random_variate drawn for the end
    of the branch).
    """
    if roeotroe*duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError("Parent rate is out of bounds, but no rate change is possible")
    r = _calc_KTB_rate(starting_rate, duration, roeotroe, rng)
    if min_rate is None:
        min_rate = 0.0
    return _bounce_constrain(starting_rate, r, min_rate, max_rate)

