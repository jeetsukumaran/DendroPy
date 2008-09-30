DendroPy Phylogenetic Computation Library

(C) 2008 Jeet Sukumaran and Mark T. Holder.

The code is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

To install the DendroPy library, end-users can run:

$ python setup.py install

Developers, who are using the "cutting-edge" code via git may prefer to use the develop command. This installs an egg-link to the path of the dendropy directory in the appropriate site-packages directory in your python installation. Thus, you will not have to reinstall after every change to the code:

$ python setup.py develop 

You might get an error message that an additional installation helper ("setuptools") needs to be downloaded:

Traceback (most recent call last):
  File "setup.py", line 29, in <module>
    from setuptools import setup, find_packages
ImportError: No module named setuptools

If so, you need to download the setuptools installer from http://peak.telecommunity.com/dist/ez_setup.py and run it using "sudo python ez_setup.py":

$ curl -O http://peak.telecommunity.com/dist/ez_setup.py
$ sudo python ez_setup.py

This will automatically download and install setuptools for you. After this you can repeat the previous step to actually install the DendroPy library:

$ sudo python setup.py install
