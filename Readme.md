NeuraPy
-------
is a collection of python modules useful for analyzing data obtained in neurophysiological experiments.

* lablib - a module for loading [lablib][lablib] data files
* monkeylogic - a module for loading data files saved by the [MonkeyLogic][ml] program
* cerebus - a module for loading .nev and .nsx files saved by BlackRock's [Cerebus][cb] software
* plexon - a module for loading .nex and .plx files produced by the Plexon [Map][map] system

[lablib]: http://maunsell.med.harvard.edu/software.html
[ml]: http://www.monkeylogic.net/
[cb]: http://www.blackrockmicro.com/content.aspx?id=13
[map]: http://www.plexon.com/product/Multichannel_Acquisition_Processor__MAP__.html

Installation
------------
1. Grab the latest code from the zip package or grab it from github by doing `git git@github.com:kghose/neurapy.git`
2. Add the directory containing `neurapy` to your modules path (e.g by adding a line to one of your .pth files)
3. Some of the modules make use of the [matplotlib][mat] package, which should be installed

[mat]: http://matplotlib.sourceforge.net/

Usage
-----

Typically, from your script, you will import one of the sub-modules e.g

`from neurapy.lablib import lablib as ll`

For module specific instructions look at the readme files in each module's directory and
look at the documentation in the source