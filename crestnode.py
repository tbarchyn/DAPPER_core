# crestnode class
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

# Purpose: this is a generic crestnode, which has advance rules and
#  internal properties. This can also be considered a 'dune slice'

# Conventions:
# 1) azimuths are in degrees
# 2) distances are in meters

import math
import ogr
from geopy import Point
from geopy.distance import VincentyDistance

from environment import *           # import environment object class

# Crestnode advance parameters
node_advance_model = 'conic'        # node advance model specification
repose = 34.0                       # angle of repose in degrees


class crestnode:
    """Crest node object: contains properties of each crestnode and advance model"""
    def __init__(self, initial_geom):
        # initial_geom is an point geometry
        self.geom = initial_geom
        self.slipface_normal = None
        self.slipface_ht = None
        self.slipface_g = None
        self.veg = None
        self.b_up = None
        
        # internal private variables for debugging
        self.q_vol = None
        self.q_az = None
        self.veg_threshold = None
        self.w = None
        self.phi = None
        self.b_down = None
        self.dh_dt_slipface = None
        self.movedist = None
        self.flux_down_slipface = None
        self.sed_delivered = None
        self.slipface_area = None
        self.move_az = None
        
        # conic slipface area variables
        self.gatept_L_lon = None
        self.gatept_L_lat = None
        self.gatept_R_lon = None
        self.gatept_R_lat = None
        self.r_inner = None
        self.r_outer = None
        self.curve_centerpoint_lon = None
        self.curve_centerpoint_lat = None
        self.delivery_gate = None
        self.open_angle = None
        
        return
        
    def advance(self, input_environment, timestep, move_for_real = True):
        # method to drive individual nodes downwind
        # this is a generic catch, where there is a choice of
        # individual node advance models
        # the move_for_real variable is a boolean, in the bluff_advance
        # method called below, the advance function is called, but the nodes
        # are not actually moved. This is so that we can update all the
        # internal variables in a fake call (and keep code simpler).
        
        # grab the parameters from the environment object
        self.q_vol = input_environment.get_sedvol(timestep)
        self.q_az = input_environment.get_dw_az(timestep)
        self.veg_threshold = input_environment.get_veg_threshold(timestep)
        
        # set the initial vegetation globally if first timestep
        if timestep == 0:
            self.veg = input_environment.get_initial_veg(timestep)
        
        # ----------------------------------------------------------------
        # Calculate the phi based on the slipface orientation and the q_az
        # calculate the angle diff between flux and the slipface normal
        if self.q_az > self.slipface_normal:
            anglespread = self.q_az - self.slipface_normal
        else:
            anglespread = self.slipface_normal - self.q_az
        
        # determine whether the angle spread is opposite to slipface direction
        # and correct anglespread to an offset angle (absolute value)
        if anglespread > 270.0:
            self.flux_down_slipface = True
            anglespread = 360.0 - anglespread
        elif anglespread <= 270.0 and anglespread > 180.0:
            self.flux_down_slipface = False
            anglespread = 360.0 - anglespread
        elif anglespread <= 180.0 and anglespread >= 90.0:
            self.flux_down_slipface = False
        elif anglespread < 90.0:
            self.flux_down_slipface = True
        else:
            print('ERROR: problem with the angle spread calc')

        # phi is the angle between the brink and the flux, 90 degrees is
        # perpendicular and should result in maximum flux
        self.phi = 90.0 - anglespread
        if self.phi < 0.0:
            self.phi = -1.0 * self.phi        # make it an absolute value

        # ----------------------------------------------------------------
        # Calculate the node advance with the node advance model
        if node_advance_model == 'trapezoid':
            self.advance_trapezoid()
        elif node_advance_model == 'conic':
            self.advance_conic()
        else:
            print('ERROR: specify node advance model in header of crestnode.py!')
        
        # ----------------------------------------------------------------
        # Calculate the movedist
        self.calc_movedist()        
        
        # Call the vegetator to check the signature of life rules
        self.vegetator()
        
        # and finally . . . move the node in space
        if move_for_real:
            self.move(self.movedist, self.move_az)
        
        return    
    
    def vegetator(self):
        # method to calculate the move distance and update the vegetation
        # Case where deposition rates exceed deposition tolerance during advance
        
        # NOTES: this is a very contentious piece of code: revist this because 
        #        this code will control the model behaviour extensively!
        
        # vegetation parameters
        max_veg_height = 2.0            # height in meters
        veg_kill_threshold = 0.5        # height in meters
        immobilization_threshold = 0.7  # height in meters where immobilization occurs
        
        if self.dh_dt_slipface > self.veg_threshold:
            if self.veg < veg_kill_threshold:
                self.veg = 0.0  # kill it all!
            else:
                # else, wear it down slowly
                veg_subtract = self.dh_dt_slipface - self.veg_threshold
                self.veg = self.veg - veg_subtract
                if self.veg < 0.0:
                    self.veg = 0.0

        # Opposite case where vegetation can survive
        else:
            veg_add = self.veg_threshold - self.dh_dt_slipface
            self.veg = self.veg + veg_add
            if self.veg > max_veg_height:
                self.veg = max_veg_height
        
        # modulate the precalculated movedist by the amount of vegetation
        self.movedist = self.movedist * (1.0 - (self.veg / immobilization_threshold))
        
        # guard against negative movedists
        if self.movedist < 0.0:
            self.movedist = 0.0
        
        return
        
    def calc_movedist(self):
        # method to set the movedist variable (need subsequent call by
        # vegetator if vegetation is incorporated in to the model)
        
        # see notebook for derivations
        
        # calc the normal movedist (down the slipface)
        normal_movedist = self.dh_dt_slipface / math.tan(math.radians(repose))
        
        # Option 1: move the crestnode in the direction of avalanches
        #self.movedist = normal_movedist
        #self.move_az = self.slipface_normal
        
        # Option 2: move the crestnode downwind, slide it along the crest
        #self.movedist = normal_movedist / math.sin(math.radians(self.phi))
        #self.move_az = self.q_az
        
        # Option 3: move the crestnode downwind based on the x component
        # of advance (see notebook for derivation)
        self.movedist = normal_movedist * math.sin(math.radians(self.phi))
        self.move_az = self.q_az
        
        # Enforce no backwards moving slipfaces?
        #if not self.flux_down_slipface:
        #    self.movedist = 0.0
            
        return
    
    def advance_trapezoid(self):
        # method to use the trapezoid method to calculate the slipface area
        # calculate b_down (bottom of trapezoid)
        self.b_down = self.b_up - 2.0 * self.w * math.tan(math.radians(self.slipface_g / 2.0))
        
        # enforce no negative b_down: this is a possibly contentious enforcement
        # that is a result of a digitizing error. In the signature of life work
        # we just omit these cells, but here we need to calculate a celerity for
        # each cell, so we have to make this enforcement.
        if self.b_down < 0.0:
            self.b_down = 0.0
        
        # calculate the deposition rate on the slipface
        self.sed_delivered = self.q_vol * self.b_up * math.sin(math.radians(self.phi))
        self.slipface_area = 0.5 * self.w * (self.b_up + self.b_down)
        self.dh_dt_slipface = self.sed_delivered / self.slipface_area

        return

    def advance_conic(self):
        # method for conical method to calculate the slipface area and
        # deposition rates. Most of the variables are pre-set within the calc_curvature
        # method of the crest object.
        if self.curve_centerpoint_lon == None or self.curve_centerpoint_lat == None:
            # the slipface is a plane, no curvature
            self.slipface_area = self.delivery_gate * self.w
        
        else:
            # else, we have a curved slipface and we can use the donut slice method to calc area
            inner_area = math.pi * self.r_inner * self.r_inner
            outer_area = math.pi * self.r_outer * self.r_outer
            donut_area = outer_area - inner_area
            donut_slice = (self.open_angle / 360.0) * donut_area
            self.slipface_area = donut_slice
                            
        self.sed_delivered = self.q_vol * self.delivery_gate * math.sin(math.radians(self.phi))
        self.dh_dt_slipface = self.sed_delivered / self.slipface_area
        
        return
    
    def move(self, distance, az):
        # a convenience method to actually move the crestnode through space
        # distance = distance to move in meters
        # az = azimuth to move the node
        lon1, lat1, z1 = self.geom.GetPoint(0)
        
        # set a geopy 'point' object and convert distance
        geopy_point1 = Point(longitude = lon1, latitude = lat1)
        distance = distance / 1000.0
        
        #TEMPORARY TEST!
        # Make the dune go backwards!
        #az = az + 180.0
        #if az > 360.0:
        #    az = az - 360.0
        
        # limit distance
        #if distance > 0.0002:
        #    distance = 0.0002
        #END TEMPORARY TEST!
        
        # compute 'vincenty' distance
        geopy_point2 = VincentyDistance(kilometers = distance).destination(geopy_point1, az)
        
        # assign new geometry to self.geom
        self.geom.SetPoint(0, x = geopy_point2.longitude, y = geopy_point2.latitude, z = z1)
        return
    
    def bluff_advance(self, input_environment, timestep):
        # a method to run the advance code without actually moving the node
        # this is used in debug operations to update all the crestnode internals
        # to be current with the positioning of the crestnode and the crest geometry
        self.advance(input_environment, timestep, move_for_real = False)
        return



