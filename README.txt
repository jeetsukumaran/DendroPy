DendroPy Phylogenetic Computation Library
=========================================

(C) 2008 Jeet Sukumaran and Mark T. Holder.

The code is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

DendroPy installs and works out-of-the-box under Python 2.5.x. It will work under Python 2.4.x if the `ElementTrees module <http://effbot.org/downloads/#elementtree>`_ is separately installed. The library itself functions perfectly fine under Python 2.6.x, but setuptools has not yet released an easyinstall egg for this version, so the automated installation of the library is not currently functioning (see below). DendroPy will not work under 2.3 without a lot of effort, and is broken under 3.0.

To install the DendroPy library, end-users can run::

    $ python setup.py install

As mentioned above, this currently only works under Python 2.5. You can also manually install the library by one of the following:

    * Add the current location of the ``dendropy`` subdirectory to your Python path environmental variable, ``$PYTHONPATH``
    
    * Copy (or symlink) the ``dendropy`` directory to the ``site-packages`` directory of your Python installation
    
You might also want to place the scripts in the ``scripts/`` subdirectory on your system path.