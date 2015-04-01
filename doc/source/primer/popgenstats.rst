*************************************
Population Genetic Summary Statistics
*************************************

The :mod:`popgenstat` module provides functions that calculate some common population genetic summary statistics.

For example, given a |DnaCharacterMatrix| as an argument, the :func:`~dendropy.popgenstat.num_segregating_sites()` function returns the raw number of segregating sites, :func:`~dendropy.popgenstat.average_number_of_pairwise_differences()` returns the average number of pairwise differences, and :func:`~dendropy.popgenstat.nucleotide_diversity()` returns the nucleotide diversity.

More complex statistics are provided by the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` class.
Objects of this class are instantatiated with two lists of |CharacterDataVector| objects as arguments, each representing a sample of DNA sequences drawn from two distinct but related populations.
Once instantiated, the following attributes of the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` object are available:

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences`
            The average number of pairwise differences between every sequence across both populations.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_between`
            The average number of pairwise differences between every sequence between both populations.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_within`
            The average number of pairwise differences between every sequence within each population.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.average_number_of_pairwise_differences_net`
            The net number of pairwise differences.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.num_segregating_sites`
            The number of segregating sites.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.wattersons_theta`
            Watterson's theta.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.wakeleys_psi`
            Wakeley's psi.

        :attr:`~dendropy.popgenstat.PopulationPairSummaryStatistics.tajimas_d`
            Tajima's D.

The following example calculates the suite of population genetic summary statistics for sequences drawn from two populations of sticklebacks.
The original data consists of 23 sequences, with individuals from Eastern Pacific populations identified by their taxon labels beginning with "``EPAC``" and individuals from Western Pacific populations identified by their taxon labels beginning with "``WPAC``".
The taxon labels thus are used as the basis for sorting the sequences into the required lists of |CharacterDataVector| objects, ``p1`` and ``p2``.

.. literalinclude:: /examples/pgstats1.py
    :linenos:

Lines 6-12 build up the two lists of |CharacterDataVector| objects by sorting the original sequences into their source populations based on the taxon label (with operational taxonomic units with labels beginning with "``EPAC``" coming from the Eastern Pacific, and assigned to the list ``p1``, while those that begin with "``WPAC``" coming from the Western Pacific, and assigned to the list ``p2``).
These lists are then passed as the instantiation arguments to the :class:`~dendropy.popgenstat.PopulationPairSummaryStatistics` constructor in line 14.
The calculations are performed immediately, and the results reported in the following lines.
