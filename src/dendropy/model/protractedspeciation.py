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
import itertools
from dendropy.utility import GLOBAL_RNG
from dendropy.utility.error import ProcessFailedException
from dendropy.utility.error import TreeSimTotalExtinctionException
from dendropy.calculate import probability

def _D(speciation_initiation_rate,
       speciation_completion_rate,
       incipient_species_extinction_rate):
    r"""
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
        species, $lambda_1$.
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
    r"""
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
        species, $lambda_1$.
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
    r"""
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
        species, $lambda_1$.
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
    r"""
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
        species, $lambda_1$.
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
    r"""
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
        species, $lambda_1$.
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
    r"""
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
        species, $lambda_1$.
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

    def __init__(self,
            speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate,
            orthospecies_extinction_rate,
            incipient_species_extinction_rate,
            lineage_label_format_template=None,
            species_label_format_template=None,
            rng=None,
            **kwargs
            ):
        self.speciation_initiation_from_orthospecies_rate = speciation_initiation_from_orthospecies_rate
        self.orthospecies_extinction_rate = orthospecies_extinction_rate
        self.speciation_initiation_from_incipient_species_rate = speciation_initiation_from_incipient_species_rate
        self.speciation_completion_rate = speciation_completion_rate
        self.incipient_species_extinction_rate = incipient_species_extinction_rate
        self.is_initial_lineage_orthospecies = kwargs.get("is_initial_lineage_orthospecies", False)
        self.species_lineage_sampling_scheme = kwargs.get("species_lineage_sampling_scheme", "random") # 'random', 'oldest', 'youngest'
        if lineage_label_format_template is None:
            self.lineage_label_format_template = "S{species_id}.L{lineage_id}"
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
        self.lineage_tree_to_species_tree_node_attr = "species_tree_node"
        self.species_tree_to_lineage_tree_node_attr = "lineage_tree_nodes"

    def generate_sample(self, **kwargs):
        """

        Samples from the Protracted Speciation Model process, returning a tuple:

            -   the lineage tree: this tree has all nodes/lineages, i.e. both
                "good" species as well as incipient species.
            -   the (ortho- or confirmed- or "good"-)species tree: the tree
                only has "good" species, i.e. with all incipient species
                pruned out.

        Each tip node on the lineage tree will have an attribute,
        ``species_tree_node``, which is a reference to the species tree tip to
        which the lineage tip corresponds.

        Each tip node on the species tree will have an attribute,
        ``lineage_tree_nodes``, which is a reference to an iterable of
        tip nodes on the lineage tree that are associated with
        this species tree node.

        Parameters
        ----------

        max_time : float or |None|
            Terminate and return results if this time is exceeded. If |None|,
            then do not terminate based on run time.
        min_extant_orthospecies : int or |None|
            Terminate and return results when at least this number of tips are
            found in the confirmed-species tree (i.e., the pruned tree
            consisting of only "good" species) and other specified conditions
            are met.
        max_extant_orthospecies : int or |None|
            Terminate and return results when at least this number of tips are
            found in the confirmed-species tree (i.e., the pruned tree
            consisting of only "good" species).
        num_extant_lineages : int or |None|
            Terminate and return results when this number of tips are found in
            the lineage tree (i.e. the tree with both incipient and good
            species). If |None|, then do not terminate based on the
            number of tipes on the incipient species tree.
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
                lineage_tree, orthospecies_tree = self._generate_trees(**kwargs)
                break
            except ProcessFailedException:
                if not is_retry_on_total_extinction:
                    raise
                num_retries += 1
                if max_retries is not None and num_retries > max_retries:
                    raise
        assert lineage_tree is not None
        return lineage_tree, orthospecies_tree

    def _generate_trees(self, **kwargs):
        max_time = kwargs.get("max_time", None)
        num_extant_lineages = kwargs.get("num_extant_lineages", None)
        min_extant_lineages = kwargs.get("min_extant_lineages", None)
        max_extant_lineages = kwargs.get("max_extant_lineages", None)
        num_extant_orthospecies = kwargs.get("num_extant_orthospecies", None)
        min_extant_orthospecies = kwargs.get("min_extant_orthospecies", None)
        max_extant_orthospecies = kwargs.get("max_extant_orthospecies", None)
        lineage_taxon_namespace = kwargs.get("lineage_taxon_namespace", None)
        species_taxon_namespace = kwargs.get("species_taxon_namespace", None)
        lineage_data = []
        for idx in range(2):
            lineage_data.append({
                "lineage_id": 0,
                "species_id": 0,
                "lineage_collection": []
                })
        phase_idx = 0
        while phase_idx < 2:
            # Run two passes if 'max_time' specified, with each pass building
            # up one side of the root and only accepted if there are at least
            # one surviving lineage; this way, condition crown age on
            # 'max_time'.
            self._generate_lineages(
                    lineage_data=lineage_data,
                    max_time=max_time,
                    num_extant_lineages=num_extant_lineages,
                    min_extant_lineages=min_extant_lineages,
                    max_extant_lineages=max_extant_lineages,
                    num_extant_orthospecies=num_extant_orthospecies,
                    min_extant_orthospecies=min_extant_orthospecies,
                    max_extant_orthospecies=max_extant_orthospecies,
                    phase_idx=phase_idx,
                    lineage_taxon_namespace=lineage_taxon_namespace,
                    species_taxon_namespace=species_taxon_namespace,
                    )
            if num_extant_orthospecies is not None or max_extant_orthospecies is not None or min_extant_orthospecies is not None:
                if "lineage_tree" in lineage_data[phase_idx]:
                    return lineage_data[phase_idx]["lineage_tree"], lineage_data[phase_idx]["orthospecies_tree"]
                else:
                    raise ProcessFailedException()
            elif num_extant_lineages is not None or max_extant_lineages is not None or min_extant_lineages is not None:
                break
            elif phase_idx == 0 and (len(lineage_data[0]["orthospecies_lineages"]) + len(lineage_data[0]["incipient_species_lineages"]) > 0):
                phase_idx += 1
            elif phase_idx == 1 and self._check_good(
                    orthospecies_lineages=lineage_data[1]["orthospecies_lineages"],
                    incipient_species_lineages=lineage_data[1]["incipient_species_lineages"]):
                phase_idx += 1
        # lineage_collection = itertools.chain(
        #         lineage_data[0]["lineage_collection"],
        #         lineage_data[1]["lineage_collection"],
        #         )
        lineage_collection = lineage_data[0]["lineage_collection"] + lineage_data[1]["lineage_collection"]
        lineage_tree = self._compile_lineage_tree(
                lineage_collection=lineage_collection,
                max_time=max_time if max_time is not None else lineage_data[phase_idx]["final_time"],
                is_drop_extinct=True,
                )
        orthospecies_tree = self._compile_species_tree(
                lineage_collection=lineage_collection,
                max_time=max_time if max_time is not None else lineage_data[phase_idx]["final_time"],
                )
        return self._finalize_trees(
                lineage_tree=lineage_tree,
                lineage_taxon_namespace=lineage_taxon_namespace,
                orthospecies_tree=orthospecies_tree,
                species_taxon_namespace=species_taxon_namespace,
                lineage_collection=lineage_collection,
                )

    def _generate_lineages(self, **kwargs):
        current_time = 0.0
        lineage_data = kwargs.get("lineage_data")
        max_time = kwargs.get("max_time", None)
        num_extant_lineages = kwargs.get("num_extant_lineages", None)
        min_extant_lineages = kwargs.get("min_extant_lineages", None)
        max_extant_lineages = kwargs.get("max_extant_lineages", None)
        num_extant_orthospecies = kwargs.get("num_extant_orthospecies", None)
        min_extant_orthospecies = kwargs.get("min_extant_orthospecies", None)
        max_extant_orthospecies = kwargs.get("max_extant_orthospecies", None)
        phase_idx = kwargs.get("phase_idx")
        lineage_taxon_namespace = kwargs.get("lineage_taxon_namespace", None)
        species_taxon_namespace = kwargs.get("species_taxon_namespace", None)
        if phase_idx == 0:
            lineage_data[phase_idx]["lineage_id"] = 1
            lineage_data[phase_idx]["species_id"] = 1
            initial_lineage = self._new_lineage(
                    lineage_id=lineage_data[phase_idx]["lineage_id"],
                    parent_lineage=None,
                    origin_time=-1e-10,
                    )
            lineage_data[phase_idx]["orthospecies_lineages"] = [initial_lineage]
            lineage_data[phase_idx]["incipient_species_lineages"] = []
        else:
            lineage_data[phase_idx]["lineage_id"] = lineage_data[0].get("lineage_id", 0) + 1
            lineage_data[phase_idx]["species_id"] = lineage_data[0].get("species_id", 0)
            initial_lineage = self._new_lineage(
                    lineage_id=lineage_data[phase_idx]["lineage_id"],
                    parent_lineage=lineage_data[0]["lineage_collection"][0],
                    origin_time=current_time,
                    )
            lineage_data[phase_idx]["orthospecies_lineages"] = []
            lineage_data[phase_idx]["incipient_species_lineages"] = [initial_lineage]
        lineage_data[phase_idx]["lineage_collection"] = [initial_lineage]
        lineage_collection = lineage_data[phase_idx]["lineage_collection"]
        orthospecies_lineages = lineage_data[phase_idx]["orthospecies_lineages"]
        incipient_species_lineages = lineage_data[phase_idx]["incipient_species_lineages"]
        while True:
            num_orthospecies = len(orthospecies_lineages)
            num_incipient_species = len(incipient_species_lineages)
            if num_incipient_species + num_orthospecies == 0:
                raise TreeSimTotalExtinctionException()
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
            # waiting_time = -math.log(self.rng.uniform(0, 1))/rate_of_any_event
            if max_time and (current_time + waiting_time) > max_time:
                current_time = max_time
                break
            # we do this here so that the (newest) tip lineages have the
            # waiting time to the next event branch lengths
            if (num_extant_lineages is not None
                    or min_extant_lineages is not None
                    or max_extant_lineages is not None):
                has_lineage_count_requirements = True
                if (
                        (num_extant_lineages is None or ((num_incipient_species + num_orthospecies) == num_extant_lineages))
                        and (min_extant_lineages is None or ((num_incipient_species + num_orthospecies) >= min_extant_lineages))
                        and (max_extant_lineages is None or ((num_incipient_species + num_orthospecies) == max_extant_lineages))
                        ):
                    is_lineage_count_requirements_met = True
                else:
                    is_lineage_count_requirements_met = False
            else:
                has_lineage_count_requirements = False
                is_lineage_count_requirements_met = None
            if max_extant_lineages is not None and (num_incipient_species + num_orthospecies) > max_extant_lineages:
                raise ProcessFailedException()
            if num_extant_orthospecies is not None or max_extant_orthospecies is not None or min_extant_orthospecies is not None:
                ## note: very expensive operation to count orthospecies leaves!
                has_orthospecies_count_requirements = True
                is_orthospecies_count_requirements_met = False
                final_time = current_time + self.rng.uniform(0, waiting_time)
                lineage_collection_snapshot = [lineage.clone() for lineage in itertools.chain(lineage_data[0]["lineage_collection"], lineage_data[1]["lineage_collection"])]
                try:
                    orthospecies_tree = self._compile_species_tree(
                            lineage_collection=lineage_collection_snapshot,
                            max_time=final_time,
                            )
                    num_leaves = len(orthospecies_tree.leaf_nodes())
                    if (
                            (num_extant_orthospecies is None or num_leaves == num_extant_orthospecies)
                            and (min_extant_orthospecies is None or num_leaves >= min_extant_orthospecies)
                            and (max_extant_orthospecies is None or num_leaves <= max_extant_orthospecies)
                            ):
                        lineage_tree = self._compile_lineage_tree(
                            lineage_collection=lineage_collection_snapshot,
                            max_time=final_time,
                            is_drop_extinct=True,
                            )
                        lineage_tree, orthospecies_tree = self._finalize_trees(
                                lineage_tree=lineage_tree,
                                lineage_taxon_namespace=lineage_taxon_namespace,
                                orthospecies_tree=orthospecies_tree,
                                species_taxon_namespace=species_taxon_namespace,
                                lineage_collection=lineage_collection_snapshot,
                                )
                        lineage_data[phase_idx]["lineage_tree"] = lineage_tree
                        lineage_data[phase_idx]["orthospecies_tree"] = orthospecies_tree
                        is_orthospecies_count_requirements_met = True
                except ProcessFailedException:
                    pass
                if max_extant_orthospecies is not None and num_leaves > max_extant_orthospecies:
                    raise ProcessFailedException
            else:
                has_orthospecies_count_requirements = False
                is_orthospecies_count_requirements_met = None
            if (
                    ( (has_lineage_count_requirements and is_lineage_count_requirements_met) and (has_orthospecies_count_requirements and is_orthospecies_count_requirements_met) )
                    or ( (has_lineage_count_requirements and is_lineage_count_requirements_met) and (not has_orthospecies_count_requirements) )
                    or ( (not has_lineage_count_requirements) and (has_orthospecies_count_requirements and is_orthospecies_count_requirements_met) )
            ):
                final_time = current_time + self.rng.uniform(0, waiting_time)
                lineage_data[phase_idx]["final_time"] = final_time
                break
            else:
                # add to current time
                current_time += waiting_time
                # Select event
                event_type_idx = probability.weighted_index_choice(weights=event_rates, rng=self.rng)
                assert (event_type_idx >= 0 and event_type_idx <= 4)
                if event_type_idx == 0:
                    # Splitting of new incipient species lineage from orthospecies lineage
                    parent_lineage = self.rng.choice(orthospecies_lineages)
                    lineage_data[phase_idx]["lineage_id"] += 1
                    new_lineage = self._new_lineage(
                            lineage_id=lineage_data[phase_idx]["lineage_id"],
                            parent_lineage=parent_lineage,
                            origin_time=current_time,
                            )
                    lineage_collection.append(new_lineage)
                    incipient_species_lineages.append(new_lineage)
                elif event_type_idx == 1:
                    # Extinction of an orthospecies lineage
                    lineage_idx = self.rng.randint(0, len(orthospecies_lineages)-1)
                    orthospecies_lineages[lineage_idx].extinction_time = current_time
                    del orthospecies_lineages[lineage_idx]
                elif event_type_idx == 2:
                    # Splitting of new incipient species lineage from incipient lineage
                    parent_lineage = self.rng.choice(incipient_species_lineages)
                    lineage_data[phase_idx]["lineage_id"] += 1
                    new_lineage = self._new_lineage(
                            lineage_id=lineage_data[phase_idx]["lineage_id"],
                            parent_lineage=parent_lineage,
                            origin_time=current_time,
                            )
                    lineage_collection.append(new_lineage)
                    incipient_species_lineages.append(new_lineage)
                elif event_type_idx == 3:
                    # Completion of speciation
                    lineage_idx = self.rng.randint(0, len(incipient_species_lineages)-1)
                    lineage = incipient_species_lineages[lineage_idx]
                    lineage.speciation_completion_time = current_time
                    lineage_data[phase_idx]["species_id"] += 1
                    lineage.species_id = lineage_data[phase_idx]["species_id"]
                    orthospecies_lineages.append(lineage)
                    del incipient_species_lineages[lineage_idx]
                elif event_type_idx == 4:
                    # Extinction of an incipient_species lineage
                    lineage_idx = self.rng.randint(0, len(incipient_species_lineages)-1)
                    incipient_species_lineages[lineage_idx].extinction_time = current_time
                    del incipient_species_lineages[lineage_idx]
                else:
                    raise Exception("Unexpected event type index: {}".format(event_type_idx))

    def _new_lineage(self,
            lineage_id,
            parent_lineage,
            origin_time=None,
            ):
        if parent_lineage is None:
            parent_lineage_id = 0
            species_id = 1
            is_parent_orthospecies = None
        else:
            parent_lineage_id = parent_lineage.lineage_id
            species_id = parent_lineage.species_id
            is_parent_orthospecies = parent_lineage.speciation_completion_time is not None
        lineage = ProtractedSpeciationProcess._Lineage(
                lineage_id=lineage_id,
                parent_lineage_id=parent_lineage_id,
                is_parent_orthospecies=is_parent_orthospecies,
                origin_time=origin_time,
                speciation_completion_time=None,
                extinction_time=None,
                species_id=species_id,
                )
        lineage._check_parent_lineage = parent_lineage # for _check_good
        return lineage

    def _correlate_lineage_and_species_trees(self, lineage_collection):
        seen_lineage_nodes = set()
        seen_species_nodes = set()
        lineage_id_to_species_maps = {}
        species_id_lineage_node_collection_map = {}
        species_id_species_node_map = {}
        for lineage in lineage_collection:
            if not lineage.lineage_node:
                assert not lineage.species_node
                continue
            if lineage.lineage_node.taxon is None:
                continue
            lineage_node = lineage.lineage_node
            species_node = lineage.species_node
            assert lineage_node not in seen_lineage_nodes
            seen_lineage_nodes.add(lineage_node)
            try:
                species_id_lineage_node_collection_map[lineage_node._species_id].add(lineage_node)
            except KeyError:
                species_id_lineage_node_collection_map[lineage_node._species_id] = set([lineage_node])
            if species_node is not None:
                assert species_node not in seen_species_nodes
                seen_species_nodes.add(lineage.species_node)
                species_id_species_node_map[species_node._species_id] = species_node
                setattr(species_node, self.species_tree_to_lineage_tree_node_attr, species_id_lineage_node_collection_map[species_node._species_id])
        for species_id, lineage_nodes in species_id_lineage_node_collection_map.items():
            for nd in lineage_nodes:
                setattr(nd, self.lineage_tree_to_species_tree_node_attr, species_id_species_node_map[species_id])

    def _finalize_trees(self,
            lineage_tree,
            lineage_taxon_namespace,
            orthospecies_tree,
            species_taxon_namespace,
            lineage_collection,
            ):
        self._build_taxa(tree=lineage_tree, taxon_namespace=lineage_taxon_namespace)
        self._build_taxa(tree=orthospecies_tree, taxon_namespace=species_taxon_namespace)
        self._correlate_lineage_and_species_trees(lineage_collection=lineage_collection)
        return lineage_tree, orthospecies_tree

    def _compile_species_tree(self,
            lineage_collection,
            max_time,
            ):
        if self.species_lineage_sampling_scheme == "oldest":
            lt = sorted(lineage_collection, key=lambda x: x.origin_time)
        elif self.species_lineage_sampling_scheme == "youngest":
            lt = sorted(lineage_collection, key=lambda x: x.origin_time, reverse=True)
        elif self.species_lineage_sampling_scheme == "random":
            lt = self.rng.sample(lineage_collection, len(lineage_collection))
        else:
            raise ValueError(self.species_lineage_sampling_scheme)
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
                )
        for k,v in to_restore_species_extinction_times.items():
            k.extinction_time = v
        return t

    def _compile_lineage_tree(self,
            lineage_collection,
            max_time,
            is_drop_extinct=True,
            ):
        return self._compile_tree(
            lineage_collection=lineage_collection,
            max_time=max_time,
            tree_type="lineage",
            is_drop_extinct=is_drop_extinct,
            )

    def _compile_tree(self,
            lineage_collection,
            max_time,
            tree_type,
            is_drop_extinct=True,
            ):
        if tree_type == "lineage":
            label_template = self.lineage_label_format_template
            node_attr = "lineage_node"
        elif tree_type == "species":
            label_template = self.species_label_format_template
            node_attr = "species_node"
        else:
            raise ValueError(tree_type)
        if len(lineage_collection) == 1:
            tree = dendropy.Tree(is_rooted=True)
            label = label_template.format(species_id=1, lineage_id=0)
            tree.seed_node._taxon_label = label
            tree.seed_node._lineage_id = 0
            tree.seed_node._species_id = 1
            tree.seed_node.edge.length = max_time
            tree.seed_node._time = max_time
            setattr(lineage_collection[0], node_attr, tree.seed_node)
            return tree
        lineage_queue = self._build_lineage_queue(
                lineage_collection=lineage_collection,
                max_time=max_time,
                is_drop_extinct=is_drop_extinct,
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
                        seed_node=seed_node,
                        is_rooted=True,
                        )
                return tree
            if parent_lineage_id == 0:
                raise ValueError

    def _build_taxa(self, tree, taxon_namespace):
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        assert len(tree.taxon_namespace) == 0
        tree.taxon_namespace = taxon_namespace
        for nd in tree.leaf_node_iter():
            nd.taxon = tree.taxon_namespace.require_taxon(label=nd._taxon_label)
            del nd._taxon_label
        return tree

    def _build_lineage_queue(self,
            lineage_collection,
            max_time,
            is_drop_extinct,
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
                node._taxon_label = label
                node._lineage_id = lineage.lineage_id
                node._species_id = lineage.species_id
                setattr(lineage, node_attr, node)
                lineageq.push_lineage(lineage=lineage, is_copy=True)
            else:
                lineageq.register_lineage_reference(lineage=lineage)
        return lineageq

    def _check_good(self, orthospecies_lineages, incipient_species_lineages):
        if orthospecies_lineages:
            return True
        if not incipient_species_lineages:
            return False
        for lineage in incipient_species_lineages:
            parent_lineage = lineage._check_parent_lineage
            origin_time = lineage.origin_time
            while parent_lineage is not None and parent_lineage.lineage_id > 1:
                if parent_lineage.speciation_completion_time is not None and parent_lineage.speciation_completion_time < origin_time:
                    return True
                origin_time = parent_lineage.origin_time
                parent_lineage = parent_lineage._check_parent_lineage
        return False
