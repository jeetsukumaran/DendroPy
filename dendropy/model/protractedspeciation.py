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
Models, modeling and model-fitting of the protracted speciation, as described in::

        Etienne, R.S., Morlon, H., and Lambert, A. 2014. Estimating the
        duration of speciation from phylogenies. Evolution 2014: 2430-2440.
        doi:10.1111/evo.12433

"""

import math
import itertools
import dendropy
from dendropy.utility import GLOBAL_RNG
from dendropy.utility.error import ProcessFailedException
from dendropy.utility.error import TreeSimTotalExtinctionException
from dendropy.calculate import probability


def _D(speciation_initiation_rate,
       speciation_completion_rate,
       incipient_species_extinction_rate):
    """
    Returns value of D, as given in eq. 5 in Etienne et al.
    (2014).

    Parameters
    ----------

    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.

    Returns
    -------
    t : float
        The duration of speciation.

    """
    D = math.sqrt(
            pow(speciation_completion_rate + speciation_initiation_rate - incipient_species_extinction_rate, 2)
            + (4.0 * speciation_completion_rate * incipient_species_extinction_rate)
        )
    return D

def _phi(speciation_initiation_rate,
       speciation_completion_rate,
       incipient_species_extinction_rate):
    """
    Returns value of $\varphi$, as given in eq. 6 in Etienne et al.
    (2014).

    Parameters
    ----------

    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.

    Returns
    -------
    t : float
        The duration of speciation.

    """
    phi = speciation_completion_rate - speciation_initiation_rate + incipient_species_extinction_rate
    return phi

def expected_duration_of_speciation(
        speciation_initiation_rate,
        speciation_completion_rate,
        incipient_species_extinction_rate,
        D=None,
        ):
    """
    Returns mean duration of speciation, following Eqs. 4 in Etienne et al.
    (2014):

        The duration of speciation differs from the speciation-completion
        time in that the latter is the waiting time until a single
        incipient lineage completes the speciation process if extinction
        was zero, whereas the former is the time needed for an incipient
        species or one of its descendants to complete speciation, condi-
        tional on the fact that speciation completes, that is, this is the
        time taken by any species that succeeded in speciating completely.

    Parameters
    ----------

    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.
    D : float
        Value of ``D`` (as given in Eq. 5 in Etienne et al. 2014). Will be
        calculated if not specified.

    Returns
    -------
    t : float
        The duration of speciation.

    """
    if D is None:
        D = _D(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    t1 = 2.0/(D - speciation_completion_rate + speciation_initiation_rate - incipient_species_extinction_rate )
    t2 = math.log(2.0/(1+((speciation_completion_rate - speciation_initiation_rate + incipient_species_extinction_rate)/D)))
    t = t1 * t2
    return t

def probability_of_duration_of_speciation(
        tau,
        speciation_initiation_rate,
        speciation_completion_rate,
        incipient_species_extinction_rate,
        D=None,
        phi=None,
        ):
    """
    Returns probability of duration of speciation, tau, following Eqs. 6
    in Etienne et al.

    Parameters
    ----------

    tau : float
        The duration of speciation.
    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.
    D : float
        Value of ``D`` (as given in Eq. 5 in Etienne et al. 2014). Will be
        calculated if not specified.
    phi : float
        Value of ``phi`` (as given in Eq. 7 in Etienne et al. 2014). Will be
        calculated if not specified.

    Returns
    -------
    p : float
        The probability of the duration of speciation, tau.

    """
    if D is None:
        D = _D(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    if phi is None:
        phi = _phi(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    n1 = 2.0 * pow(D, 2) * math.exp(-D * tau) * (D + phi)
    d1 = pow(D + phi + math.exp(-D * tau) * (D-phi), 2)
    return n1/d1

def log_probability_of_duration_of_speciation(
        tau,
        speciation_initiation_rate,
        speciation_completion_rate,
        incipient_species_extinction_rate,
        D=None,
        phi=None,
        ):
    """
    Returns probability of duration of speciation, tau, following Eqs. 6
    in Etienne et al.

    Parameters
    ----------

    tau : float
        The duration of speciation.
    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.
    D : float
        Value of ``D`` (as given in Eq. 5 in Etienne et al. 2014). Will be
        calculated if not specified.
    phi : float
        Value of ``phi`` (as given in Eq. 7 in Etienne et al. 2014). Will be
        calculated if not specified.

    Returns
    -------
    p : float
        The probability of the duration of speciation, tau.

    """
    if D is None:
        D = _D(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    if phi is None:
        phi = _phi(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    n1 = math.log(2.0) + (2 * math.log(D)) - (D * tau) + math.log(D + phi)
    d1 = 2 * (math.log(D + phi + math.exp(-D * tau)*(D-phi)))
    return n1 - d1

def maximum_probability_duration_of_speciation(
        speciation_initiation_rate,
        speciation_completion_rate,
        incipient_species_extinction_rate,
        D=None,
        phi=None,
        ):
    """
    Returns duration of speciation that maximizes probability under given
    process parameters, following eq. 8 of Etienne et al (2014).

    Parameters
    ----------

    speciation_initiation_rate : float
        The birth rate, b (the incipient species birth
        rate and the "good" species birth rate are assumed to be equal):
        the rate at which new (incipient) species are produced from
        either incipient or "good" species lineages.
    speciation_completion_rate : float
        The rate at which incipient species get converted to good
        species, $\lambda_1$.
    incipient_species_extinction_rate : float
        The incipient species exctinction rate, $\mu_1$: the rate at which
        incipient species go extinct.
    D : float
        Value of ``D`` (as given in Eq. 5 in Etienne et al. 2014). Will be
        calculated if not specified.
    phi : float
        Value of ``phi`` (as given in Eq. 7 in Etienne et al. 2014). Will be
        calculated if not specified.

    Returns
    -------
    t : float
        The duration of speciation with the maximum probability under the
        given process parameters.

    """
    if D is None:
        D = _D(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    if phi is None:
        phi = _phi(
            speciation_initiation_rate=speciation_initiation_rate,
            speciation_completion_rate=speciation_completion_rate,
            incipient_species_extinction_rate=incipient_species_extinction_rate)
    x = 1.0/D * math.log((D-phi)/(D+phi))
    return max(0, x)

class ProtractedSpeciationProcess(object):

    class ProtractedSpeciationProcessLineage(object):

        def __init__(self,
                index,
                parent_lineage,
                speciation_initiation_time,
                is_orthospecies,
                orthospecies_index):
            self._index = index
            self.is_orthospecies = is_orthospecies
            self.parent_lineage = parent_lineage
            self.speciation_initiation_time = speciation_initiation_time
            self.speciation_completion_time = None
            self.extinction_time = None
            self.lineage_tree_node_history = []
            self._label = "L{}".format(self._index)
            self.orthospecies_index = orthospecies_index

        def _get_node(self):
            return self.lineage_tree_node_history[-1]
        def _set_node(self, node):
            self.lineage_tree_node_history.append(node)
        node = property(_get_node, _set_node)

        def _get_index(self):
            return self._index
        index = property(_get_index)

        def _get_label(self):
            return self._label
        label = property(_get_label)

        def __str__(self):
            return self._get_label()

    def __init__(self,
            speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate,
            orthospecies_extinction_rate,
            incipient_species_extinction_rate,
            rng=None,
            ):
        self.speciation_initiation_from_orthospecies_rate = speciation_initiation_from_orthospecies_rate
        self.orthospecies_extinction_rate = orthospecies_extinction_rate
        self.speciation_initiation_from_incipient_species_rate = speciation_initiation_from_incipient_species_rate
        self.speciation_completion_rate = speciation_completion_rate
        self.incipient_species_extinction_rate = incipient_species_extinction_rate
        if rng is None:
            self.rng = GLOBAL_RNG
        else:
            self.rng = rng
        self.tree_factory = dendropy.Tree
        self.node_factory = dendropy.Node
        self.reset()

    def reset(self):
        self.current_time = 0.0
        self.current_lineage_index = 0
        self.current_orthospecies_index = 0
        self.current_node_index = 0
        self.current_orthospecies_lineages = []
        self.current_incipient_species_lineages = []
        self.lineage_to_orthospecies_tree_node_map = {}

    def generate_sample(self, **kwargs):
        """

        Samples from the Protracted Speciation Model process, returning a tuple:

            -   the lineage tree: this tree has all nodes/lineages, i.e. both
                "good" species as well as incipient species.
            -   the (ortho- or confirmed- or "good"-)species tree: the tree
                only has "good" species, i.e. with all incipient species
                pruned out.

        Each node on the protracted speciation tree as will as the "good" species
        tree will have an attribute, ``protracted_speciation_model_lineage``,
        which is a reference to a
        :class:`~dendropy.model.birthdeath.ProtractedSpeciationProcess.ProtractedSpeciationProcessLineage`
        instance which represents the lineage associated with this node. Note
        that each node can only be associated with a single lineage, but a
        lineage might span several nodes.

        If ``is_correlate_lineage_and_species_trees`` is |True|,
        then additional attributes will be added. See
        :meth:``dendropy.model.protractedspeciation.ProtractedSpeciationProcess.correlate_lineage_and_species_trees`
        for details.

        Parameters
        ----------

        max_time : float or |None|
            Terminate and return results when this time is reached. If |None|,
            then do not terminated based on run time.
        max_extant_orthospecies : int or |None|
            Terminate and return results when this number of tips are found in
            the confirmed-species tree (i.e., the pruned tree consisting of only
            "good" species). If |None|, then do not terminate
            based on the number of tipes on the confirmed-species tree.
        max_extant_lineages : int or |None|
            Terminate and return results when this number of tips are found in
            the lineage tree (i.e. the tree with both incipient and good
            species). If |None|, then do not terminate based on the
            number of tipes on the incipient species tree.
        is_initial_lineage_orthospecies : bool
            Whether the first lineage that initialies the process is a
            "good" species or not. Defaults to |True|: first species on
            the tree is a "good" species.
        is_retry_on_total_extinction : bool
            If |False|, then a TreeSimTotalExtinctionException will be raised
            if all lineages go extinct before the termination conditions are
            met. Defaults to |True|: if all lineages go extinct before the
            termination conditions are met, then the simulation is rerun, up to
            a maximum of ``max_retries``.
        max_retries : int
            Maximum number of runs to execute in the event of
            prematurely-terminated simulations due to all lineages going
            extinct. Once this number or re-runs is exceed, then
            TreeSimTotalExtinctionException is raised. Defaults to 1000. Set to
            |None| to never quit trying.
        is_correlate_lineage_and_species_trees: bool
            [NOT USED]
            If |True| then additional attributes will be added to the
            resulting trees to relate them. See
            :meth:``dendropy.model.protractedspeciation.ProtractedSpeciationProcess.correlate_lineage_and_species_trees`
            for details.

        Returns
        -------
        lineage_tree : |Tree| instance
            A tree from the protracted speciation process, with all lineages
            (good species as well as incipient species).
        orthospecies_tree : |Tree| instance
            A tree from the protracted speciation process with only "good" species.

        """
        is_retry_on_total_extinction = kwargs.pop("is_retry_on_total_extinction", True)
        max_retries = kwargs.pop("max_retries", 1000)
        num_retries = 0
        lineage_tree = None
        orthospecies_tree = None
        while True:
            try:
                lineage_tree, orthospecies_tree = self._run_protracted_speciation_process(**kwargs)
                break
            except ProcessFailedException:
                if not is_retry_on_total_extinction:
                    raise
                num_retries += 1
                if max_retries is not None and num_retries > max_retries:
                    raise
        assert lineage_tree is not None
        return lineage_tree, orthospecies_tree

    # def correlate_lineage_and_species_trees(self,
    #         lineage_tree,
    #         orthospecies_tree):
    #     """
    #     Correlates the protracted speciation tree and the corresponding pruned
    #     "good" species tree from a single sample of the protracted speciation
    #     process (i.e., a call to ``generate_sample()``).

    #     Each node on the lineage tree will have the following
    #     attributes added:

    #         - ``is_orthospeciation_event`` : |True| if the node represents a
    #          "good" speciation event, |False| otherwise.

    #     Each internal node on the "good" species tree will have the following
    #     attributes added:

    #         - ``lineage_tree_node``: a reference to the node on
    #           the protratcted speciation tree to which it corresponds.

    #     Each leaf node on the orthospecies tree will have the following
    #     attributes added:

    #         - ``included_lineage_tree_leaf_nodes``: the set of
    #           terminal/leaf nodes on the protracted speciation tree which are
    #           descended/included in it.

    #     """
    #     return lineage_tree, orthospecies_tree
    #     lineage_tree.calc_node_ages()

    #     orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map = {}
    #     for os_nd in orthospecies_tree.leaf_node_iter():
    #         if os_nd.parent_node is None:
    #             raise ProcessFailedException()
    #         age = str(os_nd.parent_node.age)
    #         orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map[age] = []
    #         for ch1 in os_nd.parent_node.child_node_iter():
    #             if not ch1._child_nodes:
    #                 # desc_lineages = set([sch1.protracted_speciation_model_lineage for sch1 in ch1.preorder_node_iter])
    #                 orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map[age].append((ch1, ch1.protracted_speciation_model_lineage))


    #             # desc_lineages = set([sch1.protracted_speciation_model_lineage for sch1 in ch1.preorder_node_iter])
    #             # orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map.append( (ch1, desc_lineages) )

    #     for ln_nd in lineage_tree.leaf_node_iter():
    #         cur_nd = ln_nd
    #         if cur_nd.parent_node is None:
    #             raise ProcessFailedException()
    #         focal_nd_siblings = [snd for snd in ln_nd.parent_node.child_node_iter() if snd is not ln_nd]
    #         focal_nd_parents = [pnd for pnd in ln_nd.ancestor_iter(inclusive=False)]
    #         chain = []
    #         age_found = False
    #         found = None
    #         while True:
    #             if cur_nd is None:
    #                 print("!!! {}: current is None".format(" > ".join(chain)) )
    #                 break
    #             chain.append(cur_nd.label)
    #             if cur_nd.parent_node is None:
    #                 print("!!! {}: parent is None".format(" > ".join(chain)) )
    #                 break
    #             age = str(cur_nd.parent_node.age)
    #             if age in orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map:
    #                 age_found = True
    #                 cur_nd_desc_lineages = set([sch.protracted_speciation_model_lineage for sch in cur_nd.preorder_iter()])
    #                 for orthospecies_nd, orthospecies_lineage in orthospecies_tree_leaf_node_parent_age_to_desc_node_and_lineages_map[age]:
    #                     if orthospecies_lineage in cur_nd_desc_lineages:
    #                         print(":::: {} found in desc of {}: {}".format(orthospecies_nd.protracted_speciation_model_lineage.label, ln_nd.label, [x.label for x in cur_nd_desc_lineages]))
    #                         orthospecies_nd.included_lineage_tree_leaf_nodes.add(ln_nd)
    #                         found = orthospecies_nd
    #                         break
    #                 break
    #             cur_nd = cur_nd.parent_node
    #         if not age_found:
    #             print("Age not found: {}".format(ln_nd.label))
    #         elif not found:
    #             print("Not found: {}".format(ln_nd.label))
    #         else:
    #             print("Found: {} => {}".format(ln_nd.label, found.label))
    #     # lineage_tree_node_to_orthospecies_leaf_map = {}
    #     # lineage_orthospecies_leaf_map = {}
    #     # for orthospecies_tree_nd in orthospecies_tree:
    #     #     if not orthospecies_tree_nd._child_nodes:
    #     #         lineage_orthospecies_leaf_map[orthospecies_tree_nd.protracted_speciation_model_lineage] = orthospecies_tree_nd
    #     #         continue
    #     #     lineage_tree_node = sorted(orthospecies_tree_nd.protracted_speciation_model_lineage.lineage_tree_node_history, key=lambda x: x.age)[-1]
    #     #     for history_nd in orthospecies_tree_nd.protracted_speciation_model_lineage.lineage_tree_node_history:
    #     #         if not orthospecies_tree_nd._child_nodes:
    #     #             lineage_tree_node_to_orthospecies_leaf_map[history_nd] = orthospecies_tree_nd
    #     #         # if history_nd.age is not None and history_nd.age >= lineage_tree_node.age:
    #     #         #     lineage_tree_node = history_nd
    #     #     lineage_tree_node.is_orthospeciation_event = True
    #     #     orthospecies_tree_nd.lineage_tree_node = lineage_tree_node

    #     # for lineage_tree_leaf in lineage_tree.leaf_node_iter():
    #     #     nd = lineage_tree_leaf
    #     #     print("#### {}: {}".format(lineage_tree_leaf.label, nd.protracted_speciation_model_lineage.is_orthospecies))
    #     #     while nd and nd.protracted_speciation_model_lineage.speciation_completion_time is not None:
    #     #         print(">>>> {}: {}".format(lineage_tree_leaf.label, nd.label))
    #     #         nd = nd.parent_node
    #     #     if not nd:
    #     #         print("!!! Failed: {}".format(lineage_tree_leaf.label))
    #     #     elif not nd.protracted_speciation_model_lineage.is_orthospecies:
    #     #         print("!!! Failed 2: {}".format(lineage_tree_leaf.label))
    #     #     else:
    #     #         current_lineages = [desc.protracted_speciation_model_lineage for desc in nd.preorder_iter()]
    #     #         print("set for {}: {}".format(lineage_tree_leaf.label, [x.label for x in current_lineages]))
    #     #         for lineage in current_lineages:
    #     #             if lineage in lineage_orthospecies_leaf_map:
    #     #                 lineage_orthospecies_leaf_map[lineage].included_lineage_tree_leaf_nodes.append(lineage_tree_leaf)
    #     #                 break
    #     #         else:
    #     #             print("Not found: {}".format(lineage_tree_leaf.label))


    #     # for lineage_tree_leaf in lineage_tree.leaf_node_iter():
    #     #     lineage = lineage_tree_leaf.protracted_speciation_model_lineage
    #     #     while lineage is not None and not lineage.is_orthospecies:
    #     #         lineage = lineage.parent_lineage
    #     #     if lineage is None:
    #     #         ## TODO: special case
    #     #         pass
    #     #     for orthospecies_tree_nd in orthospecies_tree:
    #     #         if orthospecies_tree_nd.protracted_speciation_model_lineage is lineage:
    #     #             if orthospecies_tree_nd._child_nodes:
    #     #                 print("Internal node {} for: {}".format(orthospecies_tree_nd.label, lineage_tree_leaf.label))
    #     #             else:
    #     #                 print("OK {} for: {}".format(orthospecies_tree_nd.label, lineage_tree_leaf.label))
    #     #             break
    #     #     else:
    #     #         print("Not found: {}".format(lineage_tree_leaf.label))



    #     # for lineage_tree_leaf in lineage_tree.leaf_node_iter():
    #     #     print("....{}".format(lineage_tree_leaf.label))
    #     #     if lineage_tree_leaf in lineage_tree_node_to_orthospecies_leaf_map:
    #     #         print("516: {}: {} => {}".format(lineage_tree_node_to_orthospecies_leaf_map[lineage_tree_leaf].label, lineage_tree_leaf.label, list(x.label for x in lineage_tree_node_to_orthospecies_leaf_map[lineage_tree_leaf].included_lineage_tree_leaf_nodes)))
    #     #         lineage_tree_node_to_orthospecies_leaf_map[lineage_tree_leaf].included_lineage_tree_leaf_nodes.append(lineage_tree_leaf)
    #     #         print("518: {}: {} => {}".format(lineage_tree_node_to_orthospecies_leaf_map[lineage_tree_leaf].label, lineage_tree_leaf.label, list(x.label for x in lineage_tree_node_to_orthospecies_leaf_map[lineage_tree_leaf].included_lineage_tree_leaf_nodes)))
    #     #     else:
    #     #         print("519")
    #     #         lineage = lineage_tree_leaf.protracted_speciation_model_lineage
    #     #         while lineage is not None and lineage not in lineage_orthospecies_leaf_map:
    #     #             lineage = lineage.parent_lineage
    #     #         if lineage is None:
    #     #             print("524")
    #     #             orthospecies_tree.seed_node.included_lineage_tree_leaf_nodes.append(lineage_tree_leaf)
    #     #             print("<<< {} >>>".format(lineage_tree_leaf.label))
    #     #         else:
    #     #             print("528")
    #     #             lineage_orthospecies_leaf_map[lineage].included_lineage_tree_leaf_nodes.append(lineage_tree_leaf)

    #         # # while lineage is not None and not lineage.is_orthospecies:
    #         # #     lineage = lineage.parent_lineage
    #         # if lineage is None:
    #         #     orthospecies_tree_nd = orthospecies_tree.seed_node
    #         # else:
    #         #     orthospecies_tree_nd = self.lineage_to_orthospecies_tree_node_map[lineage]
    #         #     # for xnd in orthospecies_tree.leaf_node_iter():
    #         #     #     if xnd.protracted_speciation_model_lineage is lineage:
    #         #     #         orthospecies_tree_nd = xnd
    #         #     #         break
    #         #     # else:
    #         #     #     orthospecies_tree_nd = orthospecies_tree.seed_node
    #         # try:
    #         #     orthospecies_tree_nd.included_lineage_tree_leaf_nodes.append(lineage)
    #         # except AttributeError:
    #         #     orthospecies_tree_nd.included_lineage_tree_leaf_nodes = [lineage]

    #         # while lineage is not None and lineage not in self.lineage_to_orthospecies_tree_node_map:
    #         #     lineage = lineage.parent_lineage
    #         # if lineage not in self.lineage_to_orthospecies_tree_node_map:
    #         #     orthospecies_tree_nd = orthospecies_tree.seed_node
    #         # else:
    #         #     orthospecies_tree_nd = self.lineage_to_orthospecies_tree_node_map[lineage]
    #         # orthospecies_tree_nd.included_lineage_tree_leaf_nodes.append(lineage)

    #     return lineage_tree, orthospecies_tree

    def _run_protracted_speciation_process(self, **kwargs):
        self.reset()
        max_time = kwargs.get("max_time", None)
        max_extant_lineages = kwargs.get("max_extant_lineages", None)
        max_extant_orthospecies = kwargs.get("max_extant_orthospecies", None)
        is_correlate_lineage_and_species_trees = kwargs.get("is_correlate_lineage_and_species_trees", False)
        taxon_namespace = kwargs.get("taxon_namespace", None)
        initial_lineage = self._new_lineage(
                parent_lineage=None,
                orthospecies_index=self.current_orthospecies_index,
                is_orthospecies=kwargs.get("is_initial_lineage_orthospecies", True),
                )
        seed_node = self._new_node(lineage=initial_lineage)
        lineage_tree = self.tree_factory( taxon_namespace=taxon_namespace, seed_node=seed_node)
        lineage_tree.is_rooted = True

        while True:

            ## Draw time to next event
            event_rates = []
            num_orthospecies = len(self.current_orthospecies_lineages)
            if max_extant_orthospecies is not None:
                ## note: expensive operation to count leaves!
                try:
                    orthospecies_tree = self._assemble_orthospecies_tree(taxon_namespace=taxon_namespace)
                    num_leaves = len(orthospecies_tree.leaf_nodes())
                    if num_leaves >= max_extant_orthospecies:
                        return self._postprocess_psm_and_orthospecies_trees(
                                orthospecies_tree=orthospecies_tree,
                                lineage_tree=lineage_tree,
                                is_correlate_lineage_and_species_trees=is_correlate_lineage_and_species_trees,
                                )
                except ProcessFailedException:
                    pass

            num_incipient_species = len(self.current_incipient_species_lineages)
            if max_extant_lineages is not None and (num_incipient_species + num_orthospecies) >= max_extant_lineages:
                break

            # Event type 0
            event_rates.append(self.speciation_initiation_from_orthospecies_rate * num_orthospecies)

            # Event type 1
            event_rates.append(self.orthospecies_extinction_rate * num_orthospecies)

            # Event type 2
            event_rates.append(self.speciation_initiation_from_incipient_species_rate * num_incipient_species)

            # Event type 3
            event_rates.append(self.speciation_completion_rate * num_incipient_species)

            # Event type 4
            event_rates.append(self.incipient_species_extinction_rate * num_incipient_species)

            # All events
            rate_of_any_event = sum(event_rates)

            # Waiting time
            waiting_time = self.rng.expovariate(rate_of_any_event)
            if max_time and (self.current_time + waiting_time) > max_time:
                t = max_time - self.current_time
                for lineage in itertools.chain(self.current_orthospecies_lineages, self.current_incipient_species_lineages):
                    lineage.node.edge.length += t
                self.current_time = max_time
                break
            self.current_time += waiting_time
            for lineage in itertools.chain(self.current_orthospecies_lineages, self.current_incipient_species_lineages):
                lineage.node.edge.length += waiting_time

            # Select event
            event_type_idx = probability.weighted_index_choice(weights=event_rates, rng=self.rng)
            assert (event_type_idx >= 0 and event_type_idx <= 4)
            # print("time {}: {}, selected = {}".format(self.current_time, event_rates, event_type_idx))

            if event_type_idx == 0:
                self._process_initiation_of_speciation_from_orthospecies(lineage_tree)
            elif event_type_idx == 1:
                self._process_orthospecies_extinction(lineage_tree)
            elif event_type_idx == 2:
                self._process_initiation_of_speciation_from_incipient_species(lineage_tree)
            elif event_type_idx == 3:
                self._process_completion_of_specation(lineage_tree)
            elif event_type_idx == 4:
                self._process_incipient_species_extinction(lineage_tree)
            else:
                raise Exception("Unexpected event type index: {}".format(event_type_idx))

            if len(self.current_orthospecies_lineages) + len(self.current_incipient_species_lineages) == 0:
                raise TreeSimTotalExtinctionException()

        orthospecies_tree = self._assemble_orthospecies_tree(taxon_namespace=taxon_namespace)
        return self._postprocess_psm_and_orthospecies_trees(
                lineage_tree=lineage_tree,
                orthospecies_tree=orthospecies_tree,
                is_correlate_lineage_and_species_trees=is_correlate_lineage_and_species_trees,
                )

    def _process_initiation_of_speciation_from_orthospecies(self, tree):
        # parent_lineage = self.rng.choice(self.current_orthospecies_lineages)
        # parent_node = parent_lineage.node
        # new_lineage = self._new_lineage(
        #         parent_lineage=parent_lineage,
        #         orthospecies_index=parent_lineage.orthospecies_index,
        #         is_orthospecies=False)
        # c1 = self._new_node(lineage=parent_lineage)
        # c2 = self._new_node(lineage=new_lineage)
        # parent_node.add_child(c1)
        # parent_node.add_child(c2)
        parent_lineage = self.rng.choice(self.current_orthospecies_lineages)
        self._process_initiation_of_speciation(parent_lineage=parent_lineage)

    def _process_initiation_of_speciation_from_incipient_species(self, tree):
        # parent_lineage = self.rng.choice(self.current_incipient_species_lineages)
        # parent_node = parent_lineage.node
        # new_lineage = self._new_lineage(
        #         parent_lineage=parent_lineage,
        #         orthospecies_index=parent_lineage.orthospecies_index,
        #         is_orthospecies=False)
        # c1 = self._new_node(lineage=parent_lineage)
        # c2 = self._new_node(lineage=new_lineage)
        # parent_node.add_child(c1)
        # parent_node.add_child(c2)
        parent_lineage = self.rng.choice(self.current_incipient_species_lineages)
        self._process_initiation_of_speciation(parent_lineage=parent_lineage)

    def _process_initiation_of_speciation(self, parent_lineage):
        parent_node = parent_lineage.node
        new_lineage = self._new_lineage(
                parent_lineage=parent_lineage,
                orthospecies_index=parent_lineage.orthospecies_index,
                is_orthospecies=False)
        c1 = self._new_node(lineage=parent_lineage)
        c2 = self._new_node(lineage=new_lineage)
        parent_node.add_child(c1)
        parent_node.add_child(c2)

        # parent_node = parent_lineage.node
        # new_lineage1 = self._new_lineage(
        #         parent_lineage=parent_lineage,
        #         orthospecies_index=parent_lineage.orthospecies_index,
        #         is_orthospecies=False,
        #         add_to_current_lineages=True)
        # new_lineage2 = self._new_lineage(
        #         parent_lineage=parent_lineage,
        #         orthospecies_index=parent_lineage.orthospecies_index,
        #         is_orthospecies=parent_lineage.is_orthospecies,
        #         add_to_current_lineages=False)
        # new_lineage2.lineage_tree_node_history = parent_lineage.lineage_tree_node_history
        # c1 = self._new_node(lineage=new_lineage1)
        # c2 = self._new_node(lineage=new_lineage2)
        # parent_node.add_child(c1)
        # parent_node.add_child(c2)

    def _process_completion_of_specation(self, tree):

        lineage = self.rng.choice(self.current_incipient_species_lineages)
        self.current_incipient_species_lineages.remove(lineage)
        self.current_orthospecies_lineages.append(lineage)
        self.current_orthospecies_index += 1
        lineage.orthospecies_index = self.current_orthospecies_index
        lineage.is_orthospecies = True
        lineage.speciation_completion_time = self.current_time

        # original_lineage = self.rng.choice(self.current_incipient_species_lineages)
        # self.current_incipient_species_lineages.remove(original_lineage)
        # self.current_orthospecies_index += 1
        # converted_lineage = self._new_lineage(
        #         parent_lineage=original_lineage,
        #         orthospecies_index=self.current_orthospecies_index,
        #         is_orthospecies=True,
        #         add_to_current_lineages=False)
        # converted_lineage.orthospecies_index = self.current_orthospecies_index
        # converted_lineage.is_orthospecies = True
        # converted_lineage.speciation_completion_time = self.current_time
        # self.current_orthospecies_lineages.append(converted_lineage)
        # converted_lineage.node = original_lineage.node

        # original_lineage = self.rng.choice(self.current_incipient_species_lineages)
        # self.current_incipient_species_lineages.remove(original_lineage)
        # self.current_orthospecies_index += 1
        # converted_lineage = self._new_lineage(
        #         parent_lineage=original_lineage,
        #         orthospecies_index=self.current_orthospecies_index,
        #         is_orthospecies=True,
        #         add_to_current_lineages=False)
        # converted_lineage.orthospecies_index = self.current_orthospecies_index
        # converted_lineage.is_orthospecies = True
        # converted_lineage.speciation_completion_time = self.current_time
        # self.current_orthospecies_lineages.append(converted_lineage)
        # c1 = self._new_node(lineage=converted_lineage)
        # original_lineage.node.add_child(c1)

        # original_lineage = self.rng.choice(self.current_incipient_species_lineages)
        # self.current_incipient_species_lineages.remove(original_lineage)
        # self.current_orthospecies_index += 1
        # converted_lineage = self._new_lineage(
        #         parent_lineage=original_lineage,
        #         orthospecies_index=self.current_orthospecies_index,
        #         is_orthospecies=True,
        #         add_to_current_lineages=False)
        # converted_lineage.orthospecies_index = self.current_orthospecies_index
        # converted_lineage.is_orthospecies = True
        # converted_lineage.speciation_completion_time = self.current_time
        # self.current_orthospecies_lineages.append(converted_lineage)
        # p2 = original_lineage.node
        # c1 = self._new_node(lineage=converted_lineage)
        # c2 = self._new_node(lineage=original_lineage)
        # p2.add_child(c1)
        # p2.add_child(c2)

        # original_lineage = self.rng.choice(self.current_incipient_species_lineages)
        # self.current_incipient_species_lineages.remove(original_lineage)
        # self.current_orthospecies_index += 1
        # converted_lineage = original_lineage
        # converted_lineage.orthospecies_index = self.current_orthospecies_index
        # converted_lineage.is_orthospecies = True
        # converted_lineage.speciation_completion_time = self.current_time
        # self.current_orthospecies_lineages.append(converted_lineage)

        # lineage = self.rng.choice(self.current_incipient_species_lineages)
        # self.current_incipient_species_lineages.remove(lineage)
        # self.current_orthospecies_index += 1
        # converted_lineage = self._new_lineage(
        #         parent_lineage=lineage.parent_lineage,
        #         orthospecies_index=self.current_orthospecies_index,
        #         is_orthospecies =True)
        # converted_lineage.speciation_completion_time = self.current_time
        # self.current_orthospecies_lineages.append(converted_lineage)
        # converted_lineage.lineage_tree_node_history = list(lineage.lineage_tree_node_history)
        # old_node = lineage.node
        # c1 = self._new_node(lineage=converted_lineage)
        # c2 = self._new_node(lineage=lineage)
        # old_node.add_child(c1)
        # old_node.add_child(c2)

    def _process_orthospecies_extinction(self, tree):
        sp = self.rng.choice(self.current_orthospecies_lineages)
        sp.extinction_time = self.current_time
        self.current_orthospecies_lineages.remove(sp)
        self._make_lineage_extinct_on_phylogeny(tree, sp.node)

    def _process_incipient_species_extinction(self, tree):
        sp = self.rng.choice(self.current_incipient_species_lineages)
        sp.extinction_time = self.current_time
        self.current_incipient_species_lineages.remove(sp)
        self._make_lineage_extinct_on_phylogeny(tree, sp.node)

    def _make_lineage_extinct_on_phylogeny(self, tree, sp):
        if len(self.current_orthospecies_lineages) == 0 and len(self.current_incipient_species_lineages) == 0:
            raise TreeSimTotalExtinctionException()
        tree.prune_subtree(sp)

    def _new_lineage(self,
            parent_lineage,
            orthospecies_index,
            is_orthospecies,
            add_to_current_lineages=True):
        self.current_lineage_index += 1
        lineage_index = self.current_lineage_index
        speciation_initiation_time = self.current_time
        new_lineage = ProtractedSpeciationProcess.ProtractedSpeciationProcessLineage(
                index=lineage_index,
                parent_lineage=parent_lineage,
                speciation_initiation_time=speciation_initiation_time,
                is_orthospecies=is_orthospecies,
                orthospecies_index=orthospecies_index)
        if add_to_current_lineages:
            if is_orthospecies:
                self.current_orthospecies_lineages.append(new_lineage)
            else:
                self.current_incipient_species_lineages.append(new_lineage)
        return new_lineage

    def _new_node(self,
            lineage,
            ):
        node = self.node_factory()
        node.edge.length = 0.0
        node.protracted_speciation_model_lineage = lineage
        node.is_orthospeciation_event = False
        self.current_node_index += 1
        node.label = "{}.n{}".format(lineage.label, self.current_node_index)
        node.index = self.current_node_index
        node.annotations.add_new(name="lineage_index", value=lineage.index)
        node.annotations.add_new(name="lineage_label", value=lineage.label)
        node.annotations.add_new(name="speciation_initiation_time", value=lineage.speciation_initiation_time)
        lineage.node = node
        return node

    def _assemble_orthospecies_tree(self, taxon_namespace=None):
        lineage_set = set(self.current_incipient_species_lineages + self.current_orthospecies_lineages)
        sorted_lineages = sorted(lineage_set,
                key = lambda x: -x.speciation_initiation_time)
        self.lineage_to_orthospecies_tree_node_map = {}
        while sorted_lineages:
            lineage = sorted_lineages.pop(0)
            lineage_set.remove(lineage)
            parent_lineage = lineage.parent_lineage
            if parent_lineage is None:
                break
            if lineage.is_orthospecies:
                orthospecies_tree_node = self._require_orthospecies_tree_node(lineage=lineage)
                # try:
                #     orthospecies_tree_node = self.lineage_to_orthospecies_tree_node_map[lineage]
                # except KeyError:
                #     orthospecies_tree_node = dendropy.Node()
                #     orthospecies_tree_node.label = "L{}".format(lineage.index)
                #     orthospecies_tree_node.protracted_speciation_model_lineage = lineage
                #     self.lineage_to_orthospecies_tree_node_map[lineage] = orthospecies_tree_node
                if lineage.is_orthospecies:
                    parent_lineage.is_orthospecies = True
                # try:
                #     orthospecies_tree_parent_node = self.lineage_to_orthospecies_tree_node_map[parent_lineage]
                # except KeyError:
                #     orthospecies_tree_parent_node = dendropy.Node()
                #     orthospecies_tree_parent_node.label = "L{}".format(parent_lineage.index)
                #     orthospecies_tree_parent_node.protracted_speciation_model_lineage = parent_lineage
                #     self.lineage_to_orthospecies_tree_node_map[parent_lineage] = orthospecies_tree_parent_node
                orthospecies_tree_parent_node = self._require_orthospecies_tree_node(lineage=parent_lineage)
                orthospecies_tree_parent_node.add_child(orthospecies_tree_node)
                if parent_lineage not in lineage_set:
                    lineage_set.add(parent_lineage)
                    sorted_lineages = sorted(lineage_set,
                            key = lambda x: -x.speciation_initiation_time)

        # identify seed node
        seed_node = None
        for nd in self.lineage_to_orthospecies_tree_node_map.values():
            if nd.parent_node is None:
                seed_node = nd
                break
        if seed_node is None:
            raise ProcessFailedException()

        # create pruned tree
        orthospecies_tree = dendropy.Tree(taxon_namespace=taxon_namespace, seed_node=seed_node)
        orthospecies_tree.is_rooted = True

        # orthospecies_tree.suppress_unifurcations()
        # orthospecies_tree.set_edge_lengths_from_node_ages()

        ## assign ages to orthospecies tree
        for nd in orthospecies_tree.postorder_node_iter():
            if nd.is_leaf():
                nd.age = 0
            else:
                nd.age = self.current_time - min(ch.protracted_speciation_model_lineage.speciation_initiation_time for ch in nd.child_node_iter())
        orthospecies_tree.set_edge_lengths_from_node_ages()
        orthospecies_tree.suppress_unifurcations()
        return orthospecies_tree

    def _postprocess_psm_and_orthospecies_trees(self,
                        lineage_tree,
                        orthospecies_tree,
                        is_correlate_lineage_and_species_trees=False,
                        ):
        # lineage_tree.suppress_unifurcations()
        for nd in lineage_tree:
            nd.annotations.add_new(name="speciation_initiation_age", value=self.current_time - nd.annotations["speciation_initiation_time"].value)
            nd.annotations.add_new(name="lineage_orthospecies_index", value=nd.protracted_speciation_model_lineage.orthospecies_index)
            nd.label = "S{}.L{}.n{}".format(nd.protracted_speciation_model_lineage.orthospecies_index, nd.protracted_speciation_model_lineage.index, nd.index)
        if is_correlate_lineage_and_species_trees:
            raise NotImplementedError()
            # return self.correlate_lineage_and_species_trees(
            #         lineage_tree=lineage_tree,
            #         orthospecies_tree=orthospecies_tree)
        else:
            return lineage_tree, orthospecies_tree

    def _require_orthospecies_tree_node(self, lineage):
        try:
            node = self.lineage_to_orthospecies_tree_node_map[lineage]
        except KeyError:
            node = self.node_factory()
            node.label = "S{}".format(lineage.orthospecies_index)

            ## this might be misleading, as, technically, multiple lineages
            ## might be associated with this node; remove?
            node.protracted_speciation_model_lineage = lineage
            node.annotations.add_new(name="lineage_index", value=lineage.index)
            node.annotations.add_new(name="lineage_label", value=lineage.label)
            node.annotations.add_new(name="lineage_orthospecies_index", value=lineage.orthospecies_index)

            self.lineage_to_orthospecies_tree_node_map[lineage] = node
        return node

    def _debug_dump_lineages(self, lineages):
        sorted_lineages = sorted(lineages,
                key = lambda x: -x.speciation_initiation_time)
        for k in sorted_lineages:
            if k.parent_lineage is None:
                pi = "NA"
                pt = "NA"
            else:
                pi = k.parent_lineage.index
                pt = k.parent_lineage.is_orthospecies
            sys.stderr.write("{:10.5f} : {:4} ({}) => {} ({})\n".format(
                    k.speciation_initiation_time,
                    k.index,
                    k.is_orthospecies,
                    pi,
                    pt))
        sys.stderr.write("\n")
