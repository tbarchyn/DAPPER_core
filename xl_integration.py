# excel integration class
# Developed at the University of Calgary
# Copyright Thomas E. Barchyn, 2014, 2015

# This file is part of DAPPER.

# DAPPER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# DAPPER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with DAPPER.  If not, see <http://www.gnu.org/licenses/>.

# Purpose: this is an object to integrate with excel via the xlwings package

# Notes about xlwings:
# This is a particularly useful, but totally flakey module that links directly
# to the excel file that is open. Here is what I've learned:
# 1) the code to call an unfrozen python file is: import file; file.run()
#    where 'run' is a function in the file. It doesn't seem to work with
#    any other type of call, e.g. execfile, or whatever
# 2) the code to call a frozen python (RunFrozenPython) needs to be edited
#    within the visual basic module to execute the name of the file directly.
#    The default tries to find something in the 'build' subdirectory.
# 3) Calling the python file from excel messes up the working directory. It
#    will find the calling file, but it will not change the working directory,
#    this creates a problem where it cannot find files that are simply
#    referred to by name, not the full pathname. This is annoying because
#    just referring to the pathname is critical to making something portable.
#    The workaround is to extract the location directory of the python
#    file, and then use that to set the working directory within the python
#    program:
#
#    exec_dir = os.path.dirname(os.path.abspath(__file__))
#    os.chdir(exec_dir)


import xlwings
import numpy

class xl_link:
    """xl link class: integrates with excel file through xlwings"""
    def __init__(self, excel_file):
        # Constructor requires the name of the excel file
        # note: this must be the full filename including path
        self.wb = xlwings.Workbook(excel_file)
        
        self.controlsheet = 'control'
        self.datasheet = 'data'
        
        # populate the basic parameters for the model
        self.start_row = int(xlwings.Range(self.controlsheet, 'B3').value)
        self.end_row = int(xlwings.Range(self.controlsheet, 'B4').value)
        self.timesteps = self.end_row - self.start_row + 1
        self.interim_prefix = str(xlwings.Range(self.controlsheet, 'B8').value)
        self.topkml = str(xlwings.Range(self.controlsheet, 'B5').value)
        self.bottomkml = str(xlwings.Range(self.controlsheet, 'B6').value)
        self.oputkml = str(xlwings.Range(self.controlsheet, 'B7').value)
        
        # get debugging boolean
        debug_cell = str(xlwings.Range(self.controlsheet, 'B9').value)
        if debug_cell == 'Y' or debug_cell == 'y':
            self.debug = True
        else:
            self.debug = False
        
        self.message('EXCEL link established, passing control to model core')
        
        return
            
    def grab_time(self):
        # Grab time from the excel file
        # construct the ranges
        column = 'A'
        start_range = column + str(self.start_row)
        end_range = column + str(self.end_row)
        
        # build an output numpy array
        outarray = numpy.zeros(self.timesteps)
        
        # grab the values from the excel file and place them in the array
        vallist = xlwings.Range(self.datasheet, start_range + ':' + end_range).value
        for i in range(0, self.timesteps):
            outarray[i] = float(vallist[i][0])
        
        return outarray
        
    def grab_winddir(self):
        # Grab winddir from the excel file
        # construct the ranges
        column = 'B'
        start_range = column + str(self.start_row)
        end_range = column + str(self.end_row)
        
        # build an output numpy array
        outarray = numpy.zeros(self.timesteps)
        
        # grab the values from the excel file and place them in the array
        vallist = xlwings.Range(self.datasheet, start_range + ':' + end_range).value
        for i in range(0, self.timesteps):
            outarray[i] = float(vallist[i][0])
        
        return outarray
        
    def grab_sedvol(self):
        # Grab sediment volume from the excel file
        # construct the ranges
        column = 'G'
        start_range = column + str(self.start_row)
        end_range = column + str(self.end_row)
        
        # build an output numpy array
        outarray = numpy.zeros(self.timesteps)
        
        # grab the values from the excel file and place them in the array
        vallist = xlwings.Range(self.datasheet, start_range + ':' + end_range).value
        for i in range(0, self.timesteps):
            outarray[i] = float(vallist[i][0])
        
        return outarray
        
    def grab_veg_threshold(self):
        # Grab vegetation survival threshold from the excel file
        # construct the ranges
        column = 'F'
        start_range = column + str(self.start_row)
        end_range = column + str(self.end_row)
        
        # build an output numpy array
        outarray = numpy.zeros(self.timesteps)
        
        # grab the values from the excel file and place them in the array
        vallist = xlwings.Range(self.datasheet, start_range + ':' + end_range).value
        for i in range(0, self.timesteps):
            outarray[i] = float(vallist[i][0])
        
        return outarray
      
    def grab_oputschedule(self):
        # Grab output schedule from the excel file
        # construct the ranges
        column = 'H'
        start_range = column + str(self.start_row)
        end_range = column + str(self.end_row)
        
        # build an output numpy array
        outarray = numpy.zeros(self.timesteps, dtype = bool)
        
        # grab the values from the excel file and place them in the array
        vallist = xlwings.Range(self.datasheet, start_range + ':' + end_range).value
        for i in range(0, self.timesteps):
            if vallist[i][0] == 'Y' or vallist[i][0] == 'y':
                outarray[i] = True
            else:
                outarray[i] = False
        
        return outarray
        
    def message(self, status):
        # method to push status messages to xl workbook
        xlwings.Range(self.controlsheet, 'B18').value = status
        return
        
