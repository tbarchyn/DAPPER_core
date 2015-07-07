# environment class
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

# Purpose: this is a generic object containing information about the
#  climate at a dune field

import numpy
from xl_integration import *

class environment:
    """Generic dunefield environment class"""
    def __init__(self, input_obj):
        # in_file can either be xl_link class, or filename
        if isinstance(input_obj, xl_link):
            # populate the requisite parameters from the excel file
            self.timesteps = input_obj.timesteps
            self.sedvol = input_obj.grab_sedvol()
            self.winddir = input_obj.grab_winddir()
            self.veg_threshold = input_obj.grab_veg_threshold()            
            self.veg_initial_conditions = 0.0
            self.time = input_obj.grab_time()
            self.outputfile = input_obj.grab_oputschedule()
        
        if isinstance(input_obj, str):
            # populate the requisite input parameters from the input file
            # the file passed to this constructor is a header file which
            # contains the simulation details. The details are contain the
            # details necessary to run the simulation.
            
            # ---------------------------------------------------------------
            # Read in the simulation details file (the header)
            simdetail_file = open(input_obj, 'r')
            eof = False
            while not eof:
                line = simdetail_file.readline()
                if line == '':
                    eof = True
                line = line.strip()
                line = line[0].split('=')
                
                # check the reads and assign parameters individually
                if line[0] == 'sim_name':
                    self.sim_name = line[1]
                
                if line[0] == 'climate_file':
                    climate_file = line[1]
                    
                if line[0] == 'interim_file_prefix':
                    self.interim_file_prefix = line[1]
                    
                if line[0] == 'input_top_kml':
                    self.input_top_kml = line[1]
                    
                if line[0] == 'input_bottom_kml':
                    self.input_bottom_kml = line[1]
                    
                if line[0] == 'final_output_file':
                    self.final_output_file = line[1]
                
                if line[0] == 'veg_initial_conditions':
                    self.veg_initial_conditions = float(line[1])
                    
                if line[0] == 'debug':
                    self.debug = line[1]
            
            simdetail_file.close()
            
            # ---------------------------------------------------------------
            # Read in the climate data
            # first read the header line
            climfile = open(climate_file, 'r')
            headerline = climfile.readline()
            headerline = headerline.split()
            headerline = headerline[0].split(',')
            headerline = numpy.array(headerline)
            climfile.close()
            
            # now read in the full climate file            
            clim_array = numpy.loadtxt(climate_file, skiprows = 1, delimter = ',', dtype = 'str')
            self.timesteps = int(clim_array.shape[0])
            
            # extract out one dimensional arrays
            time_pull = clim_array[:, headerline == 'time']
            self.time = numpy.array(time_pull, dtype = 'float')
            
            sedvol_pull = clim_array[:, headerline == 'sed_vol']
            self.sedvol = numpy.array(sedvol_pull, dtype = 'float')
            
            winddir_pull = clim_array[:, headerline == 'wind_dir']
            self.winddir = numpy.array(winddir_pull, dtype = 'float')
            
            veg_threshold_pull = clim_array[:, headerline == 'veg_threshold']
            self.veg_threshold = numpy.array(veg_threshold_pull, dtype = 'float')
            
            outputfile_pull = clim_array[:, headerline == 'output_file']
            self.outputfile = numpy.array(False, dtype = 'bool')
            for i in range(0, self.timesteps):
                if outputfile_pull[i] == 'Y' or outputfile_pull[i] == 'y':
                    self.outputfile[i] = True
                else:
                    self.outputfile[i] = False
            
        return
    
    def get_sedvol(self, timestep):
        # method to extract q at a given timestep
        q_vol = self.sedvol[timestep]     
        return float(q_vol)
        
    def get_dw_az(self, timestep):
        # method to extract azimuth at a given timestep
        # note the input was 'wind direction', which is 180 degrees
        # opposite to the direction the sediment is going to move
        wind_direction = self.winddir[timestep]
        q_direction = wind_direction + 180.0
        if q_direction > 360.0:
            q_direction = q_direction - 360.0
                
        return float(q_direction)
    
    def get_veg_threshold(self, timestep):
        # method to extract vegetation survival threshold at a given timestep
        veg_threshold = self.veg_threshold[timestep]
        return float(veg_threshold)
        
    def get_initial_veg(self, timestep):
        # method to set the initial conditions for vegetation
        return self.veg_initial_conditions
    
    def get_time(self, timestep):
        # method to extract the real time as a float
        time_float = self.time[timestep]
        return time_float
        
