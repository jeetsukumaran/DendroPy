*****
PhyML
*****

The :mod:`~dendropy.interop.phyml` module contains a wrapper for running the maximum-likelihood estimation program PhyML. You can estimate a tree given a data matrix, or substitution model parameters given a tree and a data matrix. The examples below will help you get started.

Estimating a tree is easy:

.. literalinclude:: /examples/phyml_estimate_tree.py


To estimate substitution model parameters, you need both a tree and a data matrix:

.. literalinclude:: /examples/phyml_estimate_model.py

