#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
import heapq
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

    class _Lineage(object):

        @classmethod
        def from_pbd_entry(cls, pbd_lineage_entry):
            parent_lineage_id = int(pbd_lineage_entry[1])
            if parent_lineage_id < 0:
                parent_lineage_id = abs(parent_lineage_id)
                is_parent_orthospecies = False
            else:
                is_parent_orthospecies = True
            extinction_time = float(pbd_lineage_entry[4])
            if extinction_time < 0:
                extinction_time = None
            lineage = cls(
                    lineage_id=int(pbd_lineage_entry[0]),
                    parent_lineage_id=parent_lineage_id,
                    is_parent_orthospecies=is_parent_orthospecies,
                    origin_time=float(pbd_lineage_entry[2]),
                    speciation_completion_time=float(pbd_lineage_entry[3]),
                    extinction_time=extinction_time,
                    species_id=int(pbd_lineage_entry[5]),
                    )
            return lineage

        def __init__(self, **kwargs):
            self.lineage_id = kwargs.pop("lineage_id")
            self.parent_lineage_id = kwargs.pop("parent_lineage_id")
            self.is_parent_orthospecies = kwargs.pop("is_parent_orthospecies")
            self._origin_time = kwargs.pop("origin_time")
            self.speciation_completion_time = kwargs.pop("speciation_completion_time")
            self.extinction_time = kwargs.pop("extinction_time")
            self.species_id = kwargs.pop("species_id")
            self.label = kwargs.pop("label", None)
            self.time = kwargs.pop("time", None)
            self.lineage_node = kwargs.pop("lineage_node", None)
            self.species_node = kwargs.pop("species_node", None)

        def __lt__(self, o):
            # shouldn't matter, but we actually want younger lineages first,
            # so we sort in *reverse* lineage_id
            return self.lineage_id > o.lineage_id

        def clone(self):
            return self.__class__(
                lineage_id=self.lineage_id,
                parent_lineage_id=self.parent_lineage_id,
                is_parent_orthospecies=self.is_parent_orthospecies,
                origin_time=self._origin_time,
                speciation_completion_time=self.speciation_completion_time,
                extinction_time=self.extinction_time,
                species_id=self.species_id,
                label=self.label,
                time=self.time,
                lineage_node=self.lineage_node,
                species_node=self.species_node,
                )

        def _get_is_extinct(self):
            return not self.extinction_time is None
        is_extinct = property(_get_is_extinct)

        # Non-mutable as _LineageQueue uses this as a sort key in the heap Only
        # way to modify is to pop, create new and reinsert
        def _get_origin_time(self):
            return self._origin_time
        origin_time = property(_get_origin_time)

    class _LineageQueue(object):

        class LineageQueueEmptyException(Exception):
            pass

        def __init__(self):
            self._lineages = []
            heapq.heapify(self._lineages)
            self._lineage_id_lineage_entry_map = {}
            self._lineage_id_to_original_lineage_map = {}

        def new_lineage(self, **kwargs):
            lineage = ProtractedSpeciationProcess._Lineage(**kwargs)
            self.push_lineage(lineage, is_copy=False)
            return lineage

        def push_lineage(self, lineage, is_copy=True):
            # assert lineage.lineage_id not in self._lineage_id_lineage_entry_map
            # assert lineage.lineage_id not in self._lineage_id_to_original_lineage_map
            if is_copy:
                stored_lineage = lineage.clone()
            else:
                stored_lineage = lineage
            self._lineage_id_to_original_lineage_map[lineage.lineage_id] = lineage
            self._lineage_id_lineage_entry_map[stored_lineage.lineage_id] = stored_lineage
            heapq.heappush(self._lineages, (-stored_lineage.origin_time, stored_lineage))

        def register_lineage_reference(self, lineage):
            self._lineage_id_to_original_lineage_map[lineage.lineage_id] = lineage

        def has_lineage(self, lineage):
            return lineage.lineage_id in self._lineage_id_lineage_entry_map

        def get_lineage(self, lineage_id):
            return self._lineage_id_lineage_entry_map[lineage_id]

        def get_lineage_by_index(self, idx):
            return self._lineages[0][1]

        def get_original_lineage(self, lineage_id):
            return self._lineage_id_to_original_lineage_map[lineage_id]

        def pop_youngest_lineage(self):
            # try:
            #     key, lineage = heapq.heappop(self._lineages)
            # except IndexError:
            #     raise ProtractedSpeciationProcess._LineageQueue.LineageQueueEmptyException from IndexError
            if not self._lineages:
                raise ProtractedSpeciationProcess._LineageQueue.LineageQueueEmptyException
            key, lineage = heapq.heappop(self._lineages)
            # del self._lineage_id_lineage_entry_map[lineage.lineage_id]
            # del self._lineage_id_to_original_lineage_map[lineage.lineage_id]
            return lineage

        def __len__(self):
            return len(self._lineages)

        def get_initial_lineage(self):
            return self._lineages[0][1]

    _ET_LINEAGE_ID = 0
    _ET_PARENT_ID = 1
    _ET_ORIGIN_T = 2
    _ET_SPECIATION_T = 3
    _ET_EXTINCTION_TIME = 4
    _ET_SPECIES_ID = 5
    _ETX_NODE = 6

    def __init__(self,
            speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate,
            orthospecies_extinction_rate,
            incipient_species_extinction_rate,
            lineage_label_format_template=None,
            species_label_format_template=None,
            rng=None,
            ):
        self.speciation_initiation_from_orthospecies_rate = speciation_initiation_from_orthospecies_rate
        self.orthospecies_extinction_rate = orthospecies_extinction_rate
        self.speciation_initiation_from_incipient_species_rate = speciation_initiation_from_incipient_species_rate
        self.speciation_completion_rate = speciation_completion_rate
        self.incipient_species_extinction_rate = incipient_species_extinction_rate
        if lineage_label_format_template is None:
            self.lineage_label_format_template = "S{species_id}.{lineage_id}"
        else:
            self.lineage_label_format_template = lineage_label_format_template
        if species_label_format_template is None:
            self.species_label_format_template = "S{species_id}"
        else:
            self.species_label_format_template = species_label_format_template
        if rng is None:
            self.rng = GLOBAL_RNG
        else:
            self.rng = rng
        self.reset()

    def reset(self):
        self._current_time = 0.0
        self._current_orthospecies_lineages = []
        self._current_incipient_species_lineages = []
        self._next_species_id = 1
        self._lineage_collection = []

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

    def _run_protracted_speciation_process(self, **kwargs):
        self.reset()
        max_time = kwargs.get("max_time", None)
        max_extant_lineages = kwargs.get("max_extant_lineages", None)
        max_extant_orthospecies = kwargs.get("max_extant_orthospecies", None)
        orthospecies_sampling_scheme = kwargs.get("orthospecies_sampling_scheme", "random")
        lineage_taxon_namespace = kwargs.get("lineage_taxon_namespace", None)
        species_taxon_namespace = kwargs.get("species_taxon_namespace", None)
        is_initial_lineage_orthospecies = kwargs.get("is_initial_lineage_orthospecies", True)
        initial_lineage = self._new_lineage(None)
        if is_initial_lineage_orthospecies:
            initial_lineage.speciation_completion_time = 0.0
            self._current_orthospecies_lineages.append(initial_lineage)
            self._current_incipient_species_lineages = []
        while True:
            num_orthospecies = len(self._current_orthospecies_lineages)
            num_incipient_species = len(self._current_incipient_species_lineages)
            if len(self._current_orthospecies_lineages) + len(self._current_orthospecies_lineages) == 0:
                raise TreeSimTotalExtinctionException()
            if max_extant_orthospecies is not None:
                ## note: expensive operation to count leaves!
                lineage_collection_snapshot = [lineage.clone() for lineage in self._lineage_collection]
                try:
                    orthospecies_tree = self._compile_species_tree(
                            lineage_collection=lineage_collection_snapshot,
                            max_time=self._current_time,
                            sampling_scheme=orthospecies_sampling_scheme,
                            taxon_namespace=species_taxon_namespace,
                            )
                    num_leaves = len(orthospecies_tree.leaf_nodes())
                    if num_leaves >= max_extant_orthospecies:
                        lineage_tree = self._compile_lineage_tree(
                            lineage_collection=lineage_collection_snapshot,
                            max_time=self._current_time,
                            is_drop_extinct=True,
                            taxon_namespace=lineage_taxon_namespace,
                            )
                        self._correlate_lineage_and_species_trees()
                        return lineage_tree, orthospecies_tree
                except ProcessFailedException:
                    pass
            if max_extant_lineages is not None and (num_incipient_species + num_orthospecies) >= max_extant_lineages:
                break
            ## Draw time to next event
            event_rates = []
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
            if max_time and (self._current_time + waiting_time) > max_time:
                self._current_time = max_time
                break
            self._current_time += waiting_time
            # Select event
            event_type_idx = probability.weighted_index_choice(weights=event_rates, rng=self.rng)
            assert (event_type_idx >= 0 and event_type_idx <= 4)
            if event_type_idx == 0:
                self._process_initiation_of_speciation_from_orthospecies()
            elif event_type_idx == 1:
                self._process_orthospecies_extinction()
            elif event_type_idx == 2:
                self._process_initiation_of_speciation_from_incipient_species()
            elif event_type_idx == 3:
                self._process_completion_of_speciation()
            elif event_type_idx == 4:
                self._process_incipient_species_extinction()
            else:
                raise Exception("Unexpected event type index: {}".format(event_type_idx))
        orthospecies_tree = self._compile_species_tree(
                lineage_collection=self._lineage_collection,
                max_time=self._current_time,
                sampling_scheme=orthospecies_sampling_scheme,
                taxon_namespace=species_taxon_namespace,
                )
        lineage_tree = self._compile_lineage_tree(
            lineage_collection=self._lineage_collection,
            max_time=self._current_time,
            is_drop_extinct=True,
            taxon_namespace=lineage_taxon_namespace,
            )
        self._correlate_lineage_and_species_trees()
        return lineage_tree, orthospecies_tree

    def _new_lineage(self, parent_lineage):
        if parent_lineage is None:
            parent_lineage_id = 0
            species_id = 1
            is_parent_orthospecies = None
        else:
            parent_lineage_id = parent_lineage.lineage_id
            species_id = parent_lineage.species_id
            is_parent_orthospecies = parent_lineage.speciation_completion_time is not None
        lineage = ProtractedSpeciationProcess._Lineage(
                lineage_id=len(self._lineage_collection)+1,
                parent_lineage_id=parent_lineage_id,
                is_parent_orthospecies=is_parent_orthospecies,
                origin_time=self._current_time,
                speciation_completion_time=None,
                extinction_time=None,
                species_id=species_id,
                )
        self._lineage_collection.append(lineage)
        self._current_incipient_species_lineages.append(lineage)
        return lineage

    def _correlate_lineage_and_species_trees(self):
        pass

    def _process_initiation_of_speciation_from_orthospecies(self):
        parent_lineage = self.rng.choice(self._current_orthospecies_lineages)
        self._process_initiation_of_speciation(parent_lineage=parent_lineage)

    def _process_orthospecies_extinction(self):
        lineage_idx = self.rng.randint(0, len(self._current_orthospecies_lineages)-1)
        del self._current_orthospecies_lineages[lineage_idx]

    def _process_initiation_of_speciation_from_incipient_species(self):
        parent_lineage = self.rng.choice(self._current_incipient_species_lineages)
        self._process_initiation_of_speciation(parent_lineage=parent_lineage)

    def _process_completion_of_speciation(self):
        lineage_idx = self.rng.randint(0, len(self._current_incipient_species_lineages)-1)
        lineage = self._current_incipient_species_lineages[lineage_idx]
        lineage.speciation_completion_time = self._current_time
        self._next_species_id += 1
        lineage.species_id = self._next_species_id
        self._current_orthospecies_lineages.append(lineage)
        del self._current_incipient_species_lineages[lineage_idx]

    def _process_incipient_species_extinction(self):
        lineage_idx = self.rng.randint(0, len(self._current_incipient_species_lineages)-1)
        del self._current_incipient_species_lineages[lineage_idx]

    def _process_initiation_of_speciation(self, parent_lineage):
        new_lineage = self._new_lineage(parent_lineage=parent_lineage)

    def _compile_species_tree(self,
            lineage_collection,
            max_time,
            sampling_scheme="oldest",
            taxon_namespace=None,
            ):
        if sampling_scheme == "oldest":
            lt = sorted(lineage_collection, key=lambda x: x.origin_time)
        elif sampling_scheme == "youngest":
            lt = sorted(lineage_collection, key=lambda x: -x.origin_time, reverse=True)
        elif sampling_scheme == "random":
            lt = self.rng.sample(lineage_collection, len(lineage_collection))
        else:
            raise ValueError(sampling_scheme)
        seen_species_ids = set()
        to_restore_species_extinction_times = {}
        for lineage_entry in lt:
            if lineage_entry.extinction_time is None:
                if lineage_entry.species_id not in seen_species_ids:
                    seen_species_ids.add(lineage_entry.species_id)
                else:
                    to_restore_species_extinction_times[lineage_entry] = lineage_entry.extinction_time
                    lineage_entry.extinction_time = max_time # pseudo-extinction
        t = self._compile_tree(
                lineage_collection=lt,
                max_time=max_time,
                tree_type="species",
                is_drop_extinct=True,
                taxon_namespace=taxon_namespace)
        for k,v in to_restore_species_extinction_times.items():
            k.extinction_time = v
        return t

    def _compile_lineage_tree(self,
            lineage_collection,
            max_time,
            is_drop_extinct=True,
            taxon_namespace=None,
            ):
        return self._compile_tree(
            lineage_collection=lineage_collection,
            max_time=max_time,
            tree_type="lineage",
            is_drop_extinct=is_drop_extinct,
            taxon_namespace=taxon_namespace,
            )

    def _compile_tree(self,
            lineage_collection,
            max_time,
            tree_type,
            is_drop_extinct=True,
            taxon_namespace=None,
            ):
        if tree_type == "lineage":
            label_template = self.lineage_label_format_template
            node_attr = "lineage_node"
        elif tree_type == "species":
            label_template = self.species_label_format_template
            node_attr = "species_node"
        else:
            raise ValueError(tree_type)
        node_slot_idx = ProtractedSpeciationProcess._ETX_NODE
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        if len(lineage_collection) == 1:
            tree = dendropy.Tree(taxon_namespace=taxon_namespace, is_rooted=True)
            label = label_template.format(species_id=1, lineage_id=0)
            tree.seed_node.taxon = tree.taxon_namespace.require_taxon(label=label)
            tree.seed_node.edge.length = max_time
            return tree
        lineage_queue = self._build_lineage_queue(
                lineage_collection=lineage_collection,
                max_time=max_time,
                is_drop_extinct=is_drop_extinct,
                taxon_namespace=taxon_namespace,
                node_attr=node_attr,
                label_template=label_template,
                )
        while True:
            daughter_lineage = lineage_queue.pop_youngest_lineage()
            parent_lineage_id = daughter_lineage.parent_lineage_id
            try:
                parent_lineage = lineage_queue.get_lineage(parent_lineage_id)
                start_time = daughter_lineage.origin_time
                sp_comp_time = daughter_lineage.speciation_completion_time
                daughter_node = getattr(daughter_lineage, node_attr)
                end_time = daughter_node._time
                parent_sp_comp_time = parent_lineage.speciation_completion_time
                parent_node = getattr(parent_lineage, node_attr)
                parent_end_time = parent_node._time
                ch1 = parent_node
                ch1.edge.length = parent_end_time - start_time
                ch2 = daughter_node
                ch2.edge.length = end_time - start_time
                new_node = dendropy.Node()
                new_node.add_child(ch1)
                new_node.add_child(ch2)
                new_node._time = daughter_lineage.origin_time
                setattr(parent_lineage, node_attr, new_node)
            except KeyError:
                if parent_lineage_id != 0:
                    parent_lineage = lineage_queue.get_original_lineage(parent_lineage_id) #lineage_id_lineage_entry_idx_map[parent_lineage_id]
                    parent_lineage_clone = parent_lineage.clone()
                    setattr(parent_lineage_clone, node_attr, getattr(daughter_lineage, node_attr))
                    parent_sp_comp_time2 = parent_lineage.speciation_completion_time
                    if not (parent_sp_comp_time2 is not None and (parent_sp_comp_time2 < daughter_lineage.origin_time)):
                        parent_lineage_clone.speciation_completion_time = daughter_lineage.speciation_completion_time
                    lineage_queue.push_lineage(parent_lineage_clone, is_copy=False)
            if len(lineage_queue) < 2:
                if len(lineage_queue) == 0:
                    raise ProcessFailedException
                elif len(lineage_queue) == 1:
                    initial_lineage = lineage_queue.get_initial_lineage()
                else:
                    initial_lineage = None
                    while True:
                        try:
                            initial_lineage = lineage_queue.pop_youngest_lineage()
                        except ProtractedSpeciationProcess._LineageQueue.LineageQueueEmptyException as e:
                            break
                    assert initial_lineage is not None
                seed_node = getattr(initial_lineage, node_attr)
                seed_node.edge.length = initial_lineage.origin_time
                tree = dendropy.Tree(
                        taxon_namespace=taxon_namespace,
                        seed_node=seed_node,
                        is_rooted=True,
                        )
                return tree
            if parent_lineage_id == 0:
                raise ValueError

    def _build_lineage_queue(self,
            lineage_collection,
            max_time,
            is_drop_extinct,
            taxon_namespace,
            node_attr,
            label_template,
            ):
        lineageq = ProtractedSpeciationProcess._LineageQueue()
        for lineage in lineage_collection:
            if not is_drop_extinct or not lineage.is_extinct:
                node = dendropy.Node()
                if is_drop_extinct:
                    node._time = max_time
                else:
                    node._time = lineage.extinction_time if lineage.extinction_time is not None else max_time
                label = label_template.format(species_id=lineage.species_id, lineage_id=lineage.lineage_id)
                node.taxon = taxon_namespace.require_taxon(label=label)
                setattr(lineage, node_attr, node)
                lineageq.push_lineage(lineage=lineage, is_copy=True)
            else:
                lineageq.register_lineage_reference(lineage=lineage)
        return lineageq
