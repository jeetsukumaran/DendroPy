###################################
Downloading and Installing DendroPy
###################################

DendroPy can be installed directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_ using a package manager such as |pip|_ or |setuptools|_, or alternatively the source code can be downloaded and manually installed.

Installing from the GitHub Repositories
=======================================

We recommend that you install directly from the main GitHub repository using |pip|_ (which works with an `Anaconda <https://www.anaconda.com/>`_ environment as well)::

    $ python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git
    $ python3 -m pip install git+git://github.com/jeetsukumaran/DendroPy.git

If you already have DendroPy installed, you can upgrade to the latest release version by using the "``--upgrade``" flag::

    $ python3 -m pip install --upgrade git+https://github.com/jeetsukumaran/DendroPy.git
    $ python3 -m pip install --upgrade git+git://github.com/jeetsukumaran/DendroPy.git

Note: If you do not have |pip|_ installed, you should *definitely* `install it <https://pip.pypa.io/en/latest/installing.html>`_ !

Installing From the Python Package Index
========================================

DendroPy is also "easy_installable" directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_::

    $ python3 -m pip install -U dendropy

Installing via Conda
====================

    $ conda install pip
    $ pip install dendropy

Source Download and Installation
================================

The latest release of DendroPy (|version|), can be downloaded directly from here:

    |dendropy_source_archive_url|

Once downloaded, it can be installed by running:

.. parsed-literal::

    $ tar -xvzf DendroPy-|version|.tar.gz
    $ cd DendroPy-|version|
    $ python3 setup.py install

Installing the Latest Development Version
=========================================

If you want to install from a particular branch, e.g., the latest development branch, "development-main", you can use::

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git@development-main

Or::

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git@development-main

And, to update to incorporate changes as they are added to the branch::

    $ python3 -m pip install --upgrade git+https://github.com/jeetsukumaran/DendroPy.git@development-main
    $ python3 -m pip install --upgrade git+git://github.com/jeetsukumaran/DendroPy.git@development-main

Cloning the Source Code Repository
==================================

The DendroPy source code is version-controlled using |Git|_, and the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_ can be cloned by running::

    $ git clone git://github.com/jeetsukumaran/DendroPy.git

If you plan to use this repository code as you main library code, you probably want to install DendroPy in developer mode::

    $ cd DendroPy
    $ python3 -m pip install -e .

You will, of course, need to get yourself |Git|_ for the above to work:

    - `Source <http://www.kernel.org/pub/software/scm/git/git-1.6.6.tar.gz>`_
    - `OS X binaries <http://code.google.com/p/git-osx-installer/downloads/list?can=3>`_
    - `Microsoft Windows <http://code.google.com/p/msysgit/downloads/list>`_
