#! /usr/bin/env python

############################################################################
##  treecalc.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################
"""
This module has several routines for calculating tree based statistics.
Some of these functions add "optional" attributes to the tree (such as node 
    depths) that are not kept up to date by the rest of dendropy.

In general you have to explicitly call the functions again if the tree is 
    modified in any way.
"""
import math

def add_depth_to_nodes(tree, prec=0.0000001):
    """Takes an ultrametric `tree` and adds a `depth` attribute to each internal
    node.  The `depth` is the sum of edge lengths from the node to the tips.
    
    If the lengths of different paths to the node differ by more than `prec`, 
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    """
    
    node = None    
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        if len(ch) == 0:
            node.depth = 0.0
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
    if node is None:
        raise ValueError("Empty tree encountered")

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic 
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).
    
    If the lengths of different paths to the node differ by more than `prec`, 
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.
    
    As a side effect a `depth` attribute is added to the nodes of the tree.
    
    Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
    molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
    (267). 2267-2272
    """
    # the equation is given by:
    #   T = \sum_{j=2}^n (jg_j)
    #   C = T \sqrt{\frac{1}{12(n-2)}}
    #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
    # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
    #   times between consecutive (in time, not along a branch) speciation times.
    node = None
    speciation_depths = []
    n = 0
    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        n_ch = len(ch)
        if n_ch == 0:
            node.depth = 0.0
            n += 1
        elif n_ch > 2:
            raise ValueError("Polytomy encountered")
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            last_child = ch[-1]
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > prec:
                    raise ValueError("Tree is not ultrametric")
            if n_ch == 2:
                speciation_depths.append(node.depth)
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_depths.sort(reverse=True)
    g = []
    older = speciation_depths[0]
    for age in speciation_depths[1:]:
        g.append(older - age)
        older = age
    g.append(older)
    if not g:
        raise ValueError("No internal nodes found (other than the root)")
    assert(len(g) == (n - 1))
    T = 0.0
    accum = 0.0
    for i in xrange(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C



def _calc_TKP_rate(starting_rate, duration, roeotroe, rng):
    """Returns a simulated rate for the head node of a tree when:
        * the tail node has rate `starting_rate`
        * the time duration of the edge is `duration`
        * the rate of evolution of the rate of evolution is `roeotroe` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    `rng` is a random number generator.
    
    The model used to generate the rate is the one described by Thorne, Kishino,
    and Painter 1998.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of `starting_rate`
    The variance of the normal distribution for the logarithm of the ending rate 
        is the product of `duration` and `roeotroe`
    """
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        mu = math.log(starting_rate)
        return rng.lognormvariate(mu, math.sqrt(rate_var))
    return starting_rate
    
def _calc_KTB_rate(starting_rate, duration, roeotroe, rng):
    """Returns a simulated rate for the head node of a tree when:
        * the tail node has rate `starting_rate`
        * the time duration of the edge is `duration`
        * the rate of evolution of the rate of evolution is `roeotroe` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    `rng` is a random number generator.
    
    The model used to generate the rate is the one described by Kishino, Thorne,
    and Bruno 2001.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of `starting_rate`
    The variance of the normal distribution for the logarithm of the ending rate 
        is the product of `duration` and `roeotroe`
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
    return bounce_constrain(starting_rate, r, min_rate, max_rate)

def bounce_constrain(start_x, x, min_x=None, max_x=None):
    """Returns the value of variable and its mean value over a path.
    We assume that some variable started at `start_x` and moved toward `x`, but 
    has to bounce of barriers specified by `min_x` and `max_x`.
    
    `x` determines the direction and magnitude of the change.
    
    `start_x` must fall in the legal range (between the min and max). If
    `x` is also legal, then (x, (x + start_x)/2.0) will be returned reflecting
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



def simulate_continuous(node, rng, **kwargs):
    """Takes a node and a random number generator object, `rng` This function
    "evolves" a set of rates on the subtree descending from the  `node`.
    
    kwargs keys that are used are:
        `roeotroe` --  the rate of evolution of the rate of evolution. This 
            parameter that controls the degree of deviations away from the 
            molecular clock.  
        `min_rate` is the minimum rate (default None)
        `max_rate' is the maximum rate (default None), 
        `model` is a string specifying the name of the model. Currently only the
            KTB (Kishino, Thorne, Bruno) is supported
        `time_attr` is a string that specifies the name of the attribute 
            that returns the branch length in terms of time for a node. The 
            default is "edge_length"
        `val_attr` is the string that specifies the name of the attribute 
            used to hold the value that is evolving along the nodes.  The root 
            of the subtree is assumed to have this field on calling of the 
            function.  On success all nodes in the subtree will have the 
            attribute.  The default is "mutation_rate"
        `mean_val_attr` if specified this is string that gives the name of 
            attribute in each node that is mean value for the branch (default is
            None). This is filled in after time_attr and val_attr are read,
            so it is permissible to have this attribute match one of thos 
            strings (although it will make the model odd if the mean_val_attr
            is the same as the val_attr)
         `constrain_rate_mode` controls the behavior when the minimum or maximum 
            rate is simulated. The choices are "crop", and "linear_bounce"
            "crop" means that the rate is set to the most extreme value allowed.
            "linear_bounce" considers the path of evolution of rate to be a
                simple line from the ancestor's rate to the proposed rate. The
                point at which the path crosses the extreme value is determined
                and the rate is "reflected off" the limiting rate at that point.
                This causes the rate to avoid the extreme values more than a 
                simulation of small time slices that simply rejects illegal 
                rates.
    
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
    nd_iter = node.preorder_iter(node)
    # skip the first node -- it should already have a rate
    nd_iter.next()
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
        rate_func = _calc_KTB_rates_crop
    else:
        rate_func = _calc_KTB_rates_linear_bounce
    for nd in nd_iter:
        starting_rate = getattr(nd.parent_node, val_attr)
        duration = getattr(nd, time_attr)
        r, mr  = rate_func(starting_rate, duration, roeotroe, rng, min_rate, max_rate)
        setattr(nd, val_attr, r)
        if mean_val_attr:
            setattr(nd, mean_val_attr, mr)
        
