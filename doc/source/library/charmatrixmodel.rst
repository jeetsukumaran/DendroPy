***********************************************************************************
Character State Sequences and Matrices -- :mod:`dendropy.datamodel.charmatrixmodel`
***********************************************************************************

.. module:: dendropy.datamodel.charmatrixmodel

.. toctree::
    :maxdepth: 2

The :class:`CharacterMatrix` Class
==================================

A range of methods exist for importing data from another
:class:`CharacterMatrix` instance.  These vary depending on how "new" and
"existing" are treated.  A "new" sequence is a sequence in the other matrix
associated with a :class:`Taxon` object for which there is no sequence defined
in the current matrix.  An "existing" sequence is a sequene in the other
matrix associated with a :class:`Taxon` object for which there *is* a sequence
defined in the current matrix.

+---------------------------------+---------------------------------------------+--------------------------------------------+
|                                 | New Sequences: IGNORED                      | New Sequences: ADDED                       |
+=================================+=============================================+============================================+
| Existing Sequences: IGNORED     | [NO-OP]                                     | :meth:`CharacterMatrix.add_sequences()`    |
+---------------------------------+---------------------------------------------+--------------------------------------------+
| Existing Sequences: OVERWRITTEN | :meth:`CharacterMatrix.replace_sequences()` | :meth:`CharacterMatrix.update_sequences()` |
+---------------------------------+---------------------------------------------+--------------------------------------------+
| Existing Sequences: EXTENDED    | :meth:`CharacterMatrix.extend_sequences()`  | :meth:`CharacterMatrix.extend_matrix()`    |
+---------------------------------+---------------------------------------------+--------------------------------------------+

.. autoclass:: dendropy.datamodel.charmatrixmodel.CharacterMatrix
    :members:
