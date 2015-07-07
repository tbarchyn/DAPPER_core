# dunefield class
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

# Purpose: this is a generic object containing a series of individual
#  crestlines: a 'dunefield'.

import ogr
import sys
import time

from crest import *                       # import crest class
from environment import *                 # import environment object class

class dunefield:
    """Dunefield class: a number of crests in a collection referred to as a dunefield"""
    def __init__(self, in_topkml, in_bottomkml):
        # in_topkml file containing top of slipface crests
        # in_bottomkml file containing bottom of slipface crests
        self.topkml = in_topkml
        self.bottomkml = in_bottomkml

        #-----------------------------------------------------------
        # Open up the top kml file
        self.input_drv_top = ogr.GetDriverByName('KML')
        self.input_data_top = self.input_drv_top.Open(self.topkml)
        if self.input_data_top == None:
            print('ERROR: cannot open top crestline file')
            time.sleep(3)
            sys.exit(1)
            
        self.input_layer_top = self.input_data_top.GetLayer()
        self.num_crests = self.input_layer_top.__len__()
        
        # Open up the bottom kml file
        self.input_drv_bottom = ogr.GetDriverByName('KML')
        self.input_data_bottom = self.input_drv_bottom.Open(self.bottomkml)
        if self.input_data_bottom == None:
            print('ERROR: cannot open bottom crestline file')
            time.sleep(3)
            sys.exit(1)
            
        self.input_layer_bottom = self.input_data_bottom.GetLayer()

        # pull out each individual top crest into crest objects
        # and try to find the matching crest bottom. Note this
        # will fail if there is any screwups in the kml files
        # TO DO: make this a bit more robust and throw some helpful errors
        # if things aren't working correctly
        self.crests = []
        for i in range(0, self.num_crests):
            # pull out a top feature and check out its id
            top_feature = self.input_layer_top.GetFeature(i)
            crestid = top_feature.GetFieldAsInteger(0)
            
            # find the corresponding bottom crest
            for j in range(0, self.num_crests):
                bot_feature = self.input_layer_bottom.GetFeature(j)
                test_crestid = bot_feature.GetFieldAsInteger(0)
                if test_crestid == crestid:
                    break
                    
            # breaking at this point will leave the bot_feature open
            # as the feature to pass to the crest object
            newcrest = crest(top_feature, bot_feature)
            self.crests.append(newcrest)
        
        return

    def advance(self, input_environment, timestep):
        # this method advances the crests based on the node advance model
        # input_environment object is the present environment
        # timestep is the timestep for reading from the input_environment object
        
        # loop through the crests, and pass instruction down to crests
        for i in range(0, self.num_crests):
            self.crests[i].advance(input_environment, timestep)
               
        return
                
    def output(self, output_kml):
        # output the present status of the crests into a kml file
        # output_kml is the output kml file name
        out_drv = ogr.GetDriverByName('KML')
        out_data = out_drv.CreateDataSource(output_kml)
        out_layer = out_data.CreateLayer(name = output_kml, 
                                        srs = self.input_layer_top.GetSpatialRef(),
                                        geom_type = ogr.wkbLineString)

        # loop through crests and add to output datasource
        for i in range(0, self.num_crests):
            out_feature = self.crests[i].feat
            out_layer.CreateFeature(out_feature)
       
        out_data.Destroy()         # flush out file to disk
        
        print ('Output file written ' + output_kml)
        
        return
        
    def debug(self, input_environment, timestep):
        # method to pass debug command and timestep down to crests
        for i in range(0, self.num_crests):
            self.crests[i].debug(input_environment, timestep)
            
        return    
        
