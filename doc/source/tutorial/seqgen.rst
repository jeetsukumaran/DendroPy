*******
Seq-Gen
*******

|DendroPy| includes native infrastructure for phylogenetic sequence simulation
on DendroPy trees under the HKY model. Being pure-Python, however, it is a
little slow. If |SeqGen| is installed on your system, though, you can take
advantage of the :class:`dendropy.interop.seqgen.SeqGen` wrapper to efficiently
simulate sequences under a wide variety of models.  The following examples
should be enough to get started, and the class is simple and straightforward
enough so that all options should be pretty much self-documented.

.. literalinclude:: /examples/seqgen.py
