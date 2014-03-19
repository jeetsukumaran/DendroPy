###################################
Downloading and Installing DendroPy
###################################

DendroPy can be installed directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_ using a package manager such as |pip|_ or |setuptools|_, or alternatively the source code can be downloaded and manually installed.

Installing From the Python Package Index
========================================

DendroPy is "easy_installable" directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_.
If you have |pip|_ set up on your system, you can install the latest release of DendroPy by running::

    $ sudo pip install dendropy

Alternatively, if you have |setuptools|_ installed, you can run::

    $ sudo easy_install -U dendropy

Source Download and Installation
================================

The latest release of DendroPy (|version|), can be downloaded directly from here:

    |dendropy_source_archive_url|

Once downloaded, it can be installed by running:

.. parsed-literal::

    $ tar -xvzf DendroPy-|version|.tar.gz
    $ cd DendroPy-|version|
    $ sudo python setup.py install

Cloning the Source Code Repository
==================================

The DendroPy source code is version-controlled using |Git|_, and the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_ can be cloned by running::

    $ git clone git://github.com/jeetsukumaran/DendroPy.git

If you plan to use this repository code as you main library code, you probably want to install DendroPy in developer mode::

    $ cd DendroPy
    $ sudo python setup.py develop

You will, of course, need to get yourself |Git|_ for the above to work:

    - `Source <http://www.kernel.org/pub/software/scm/git/git-1.6.6.tar.gz>`_
    - `OS X binaries <http://code.google.com/p/git-osx-installer/downloads/list?can=3>`_
    - `Microsoft Windows <http://code.google.com/p/msysgit/downloads/list>`_
