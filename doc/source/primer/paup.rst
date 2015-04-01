****
PAUP
****

The :mod:`~dendropy.interop.paup` module provides functions to estimate a tree given a data matrix, or a substitution model given a tree and a data model.

Trees can be estimated using likelihood:

.. literalinclude:: /examples/paup_estimate_tree_ml.py

Or neighbor-joining:

.. literalinclude:: /examples/paup_estimate_tree_nj.py

Estimating a substitution model parameters requires both a tree and a data matrix:

.. literalinclude:: /examples/paup_estimate_model.py
