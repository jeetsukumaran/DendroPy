*****
RAxML
*****

The :mod:`~dendropy.interop.raxml` module provides the :class:`~dendropy.interop.raxml.RaxmlRunner` class, which is a lighweight (i.e., mostly "pass-through") wrapper to the |RAxML| maximum-likelihood tree estimation program. |RAxML| needs to be installed in the system for this class to be used.

The class handles the exporting of the |DendroPy| dataset in a format suitable for |RAxML| analysis, and re-reading the resulting tree back into a |DendroPy| object.

The basic call assumes nucleotide data and estimates a tree under the ``GTRCAT`` model:


.. literalinclude:: /examples/raxml_estimate_tree.py


Currently, the only way to customize the call to the underlying |RAxML| program using :meth:`~dendropy.interop.raxml.RaxmlRunner.estimate_tree` is to pass it a list of command-line arguments, with each argument token a separate list element. So, for example, to include invariant sites in the substitutio model and run 250 independent tree searches::


    >>> tree = rx.estimate_tree(data, ['-m', 'GTRCATI', '-N', '250'])

Obviously, while this works, it is neither ideal nor very Pythonic. Future releases will polish up the interface.
