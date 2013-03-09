NeuraPy
-------
is a collection of python modules useful for analyzing data obtained in neurophysiological experiments.

* cerebus - a module for loading .nev and .nsx files saved by BlackRock's [Cerebus][cb] software
* lablib - a module for loading [lablib][lablib] data files
* monkeylogic - a module for loading data files saved by the [MonkeyLogic][ml] program
* neuralynx - a module for reading files produced by [neuralynx][nl] boxes
* neuroexplorer - reads/writes .nex files for the neuroexplorer/offline sorter programs

[lablib]: http://maunsell.med.harvard.edu/software.html
[ml]: http://www.monkeylogic.net/
[cb]: http://www.blackrockmicro.com/content.aspx?id=13
[map]: http://www.plexon.com/product/Multichannel_Acquisition_Processor__MAP__.html
[nl]: http://neuralynx.com/products/digital_data_acquisition_systems/

Installation
------------
1. Grab the latest code from the zip package or grab it from github by doing `git git://github.com/kghose/neurapy.git`
2. Add the directory containing `neurapy` to your modules path (e.g by adding a line to one of your .pth files)
3. Some of the modules make use of the [matplotlib][mat] package, which should be installed

[mat]: http://matplotlib.sourceforge.net/

Usage
-----

Typically, from your script, you will import one of the sub-modules e.g

`from neurapy.lablib import lablib as ll`

For module specific instructions look at the readme files in each module's directory and
look at the documentation in the source

Licence
-------
neurapy is available for download for neurophysiologists in the hope that the programs may of use to other researchers.
No claim is made that the code works and is free of bugs, but I would appreciate bug reports (kaushik.ghose@gmail.com)

http://www.gnu.org/licenses/gpl.html

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
