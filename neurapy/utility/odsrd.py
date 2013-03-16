"""A simple reader for .ods files. Requires the great odfpy library

from neurapy.utility import odsrd
doc = odsrd.ODSReader("../../Notes/sessions_and_neurons.ods")
sh = doc.sheet_by_name('Neurons')
print sh

Parts of the code were inspired by Marco Conti's use of odfpy.
(http://www.marco83.com/work/173/read-an-ods-file-with-python-and-odfpy/)
"""

import pylab, odf.opendocument
from odf.table import Table, TableRow, TableCell
from odf.text import P

class ODSReader:

  def __init__(self, file):
    """Read metadata from file."""
    self.doc = odf.opendocument.load(file)
    self.SHEETS = {}
    self.__TABLEELEMENT__ = {}
    for sheet in self.doc.spreadsheet.getElementsByType(Table):
      self.prepareSheet(sheet)

  def prepareSheet(self, sheet):
    """Read the sheet name and put it in the list."""
    name = sheet.getAttribute("name")
    self.SHEETS[name] = None
    self.__TABLEELEMENT__[name] = sheet

  def sheet_by_name(self, name):
    """Return a sheet to us. Load on demand. Load if not already loaded"""
    if self.SHEETS[name] is None:
      self.readSheet(self.__TABLEELEMENT__[name])
    return self.SHEETS[name]


  def readSheet(self, sheet):
    """On demand, actually read a sheet and return it to us."""
    name = sheet.getAttribute("name")
    rows = sheet.getElementsByType(TableRow)
    arrRows = []
    join = ''.join

    for row in rows[:-2]:#Last two 'rows' are bogus
      arrCells = []
      cells = row.getElementsByType(TableCell)

      for cell in cells:
        ps = cell.getElementsByType(P)
        #http://www.skymind.com/~ocrow/python_string/
        #http://stackoverflow.com/questions/3766711/python-advanced-nested-list-comprehension-syntax
        cell_text = join([unicode(n.data) for p in ps for n in p.childNodes if n.nodeType==3])

        crpt = cell.getAttribute("numbercolumnsrepeated")
        if(not crpt):
          arrCells.append(cell_text)
        else:
          arrCells += [cell_text] * int(crpt)

      rrpt = row.getAttribute("numberrowsrepeated")
      if not rrpt:
        arrRows.append(arrCells)
      else:
        arrRows += [arrCells] * int(rrpt)

    self.SHEETS[name] = pylab.array(arrRows, dtype='str') #This allows us to slice the table efficiently and do finds on it

  # returns a sheet as an array (rows) of arrays (columns)
  def getSheet(self, name):
    return self.SHEETS[name]

