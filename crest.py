# crest class
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

# Purpose: this is a generic crest class for use in monitoring based
#  dune advance studies. The crestline is iteratively advanced downwind
#  with sediment transport inputs, can operate with or without vegetation

import ogr
import os
from geopy import Point
from geopy.distance import VincentyDistance

from crestnode import *                # the crestnode class
from environment import *              # import environment object class

# Crest internal parameters
repose = 34.0                          # angle of repose in degrees
extra_detailed_debugging = False       # only set to true if running 1 crest, 1 timestep


class crest:
    """Generic crest class"""
    def __init__(self, top_feature, bottom_feature):
        # constructor requires two features for the bottom and top of
        # crests. The bottom feature is just used presently to calculate
        # the height of the slipface, and then is discarded.
        self.feat = top_feature
        self.bot_feat = bottom_feature
        
        # the 'name' as set by google earth is the first field, set internally
        self.id = self.feat.GetFieldAsInteger(0)          # id tag
    
        # decompose the top of the crest into a series of crestnodes
        self.crestnodes = []
        self.geom = self.feat.GetGeometryRef()
        self.num_crestnodes = self.geom.GetPointCount()
        
        for i in range(0, self.num_crestnodes):
            new_lon, new_lat, new_z = self.geom.GetPoint(i) # as list
            new_point = ogr.Geometry(ogr.wkbPoint)
            new_point.AddPoint(x = new_lon, y = new_lat, z = new_z)
            new_crestnode = crestnode(new_point)
            self.crestnodes.append(new_crestnode)
        
        # pull out the bottom geometry
        self.bot_geom = self.bot_feat.GetGeometryRef()
        
        self.downhill_direction = None    # flag for storing downhill dir
        
        self.update_geom()                # update the geometry for the first time
        return
        
    def calc_normals_heights(self):
        # method to calculate the normals to crestnodes and evaluate heights
        # this needs to be done together because the normals need to be assigned
        # based on intersection with the bottom of slipface line.
        
        # determine angles for crestnodes by looping across each crestnode
        for i in range(0, self.num_crestnodes):
            # Run special conditions for first and last nodes
            if i == 0:
                # FIRST NODE
                # az_cn_N = azimuth to next crestnode
                az_cn_N = self.calc_azimuth(i, i+1)
                
                # set as 90 degrees to the azimuth to start
                normal_az = az_cn_N - 90.0
                if normal_az < 0.0:
                    normal_az = normal_az + 360.0
                
            elif i == (self.num_crestnodes - 1):
                # LAST NODE
                # az_cn_P = azimuth to previous crestnode
                az_cn_P = self.calc_azimuth(i, i-1)
                
                # set as 90 degrees to the azimuth to start
                normal_az = az_cn_P - 90.0
                if normal_az < 0.0:
                    normal_az = normal_az + 360.0
                    
            else:
                # REGULAR MID-CREST CALCULATIONS
                # az_cn_P = azimuth to previous crestnode
                # az_cn_N = azimuth to next crestnode
                az_cn_P = self.calc_azimuth(i, i-1)
                az_cn_N = self.calc_azimuth(i, i+1)
                
                # determine azimuth spread and determine average
                az_spread = az_cn_N - az_cn_P
                if az_spread < 0.0:
                    az_spread = -1.0 * az_spread
                
                if az_cn_N < az_cn_P:
                    normal_az = az_cn_N + (az_spread / 2.0)
                else:
                    normal_az = az_cn_P + (az_spread / 2.0)
                
                if normal_az > 360.0:
                    normal_az = normal_az - 360.0

            # -----------------------------------------------------------
            # Now we have the azimuth normal value based on the angle split
            # between the two adjacent nodes, but we aren't sure which
            # way is down yet, it could be 180 degrees out. In the initial setup
            # we will use the intersection with the bottom slipface line to
            # determine which way is down. For subsequent calls to this method
            # we will evaluate the side of the crestline which contains the
            # 'down direction', and use this to ensure continuity.
            
            # INITIAL SETUP CALL: assign slipface normal and slipface height
            if self.crestnodes[i].slipface_normal == None:
                # Here we don't know the direction which corresponds to down
                # and we must determine the direction with a ray and intersect test.
                slip_ht = self.calc_slipface_height(i, normal_az)
                if slip_ht == None:
                    normal_az = self.flip_az(normal_az)
                    slip_ht = self.calc_slipface_height(i, normal_az)
                    if slip_ht == None:
                        print('ERROR: bottom intersect problem, crestid: ',
                              str(self.id), ', crestnode: ', str(i))
                
                # OK, we should have a slipface height and set the
                # slipface normal properly to be downhill. We can now
                # populate the crestnode objects.
                self.crestnodes[i].slipface_normal = normal_az
                self.crestnodes[i].slipface_ht = slip_ht
                self.crestnodes[i].w = slip_ht / math.tan(math.radians(repose))
            
            # SUBSEQUENT CALLS: only assign slipface normal
            else:
                # Here we know which direction is downhill on the crest, based
                # on the function we called previously to set the downhill
                # direction. The purpose of this function is to check the normal
                # azimuth to set it properly based on the downhill direction of the
                # crest.
                # assign the slipface normal to start the test
                self.crestnodes[i].slipface_normal = normal_az
                if self.calc_downhill_direction(i) != self.downhill_direction:
                    # if the calculated slipface normal direction doesn't match the
                    # pre-calculated downhill direction, we have got the normal
                    # az wrong, and we need to flip it 180 degrees.
                    normal_az = self.flip_az(normal_az)
                    self.crestnodes[i].slipface_normal = normal_az
                    
                    # let's do a double check to ensure there are no issues
                    # with the functions (perhaps remove this to speed up things later).
                    if self.calc_downhill_direction(i) != self.downhill_direction:
                        print('ERROR: slipface normal assignment error!')

        return
    
    def set_downhill_direction(self):
        # method to set the internal downhill_direction flag
        # this method gets called once in setup, triggered by None 
        # set by the constructor.
        
        if self.downhill_direction == None:
            # calculate for crestnode 1 for now, perhaps expand to check
            self.downhill_direction = self.calc_downhill_direction(1)
        
        return
            
    def calc_downhill_direction(self, target_crestnode_index):           
        # method to calculate the downhill direction based on normal azimuth
        # target crestnode index can be 0 or self.num_crestnodes, the function
        # has conditional splits for these situations.
            
        # The returned downhill direction is either 'left' or 'right', which denotes
        # which side of the crest is downhill when traveling down the
        # crest in the direction of increasing crestnode index.
        
        # operate program when we have a backwards and forwards azimuths
        if target_crestnode_index > 0 and target_crestnode_index < (self.num_crestnodes - 1):
            # terminology: crestnode 1 is target, crestnode 0 is one previous, and
            # crestnode 2 is one ahead
            
            # get azimuth from crestnode 1 to crestnode 0
            az_1_0 = self.calc_azimuth(target_crestnode_index, (target_crestnode_index - 1))
            
            # get azimuth from crestnode 1 to crestnode 2
            az_1_2 = self.calc_azimuth(target_crestnode_index, (target_crestnode_index + 1))
            
            # get slipface normal from the crestnode (which should be correct)
            slipface_normal = self.crestnodes[target_crestnode_index].slipface_normal
            
            # set conditions based on relative ordering of azimuths
            if az_1_2 > slipface_normal and slipface_normal > az_1_0:
                direction = 'left'
            elif slipface_normal > az_1_2 and az_1_2 > az_1_0:
                direction = 'right'
            elif az_1_2 > az_1_0 and az_1_0 > slipface_normal:
                direction = 'right'
            elif az_1_0 > az_1_2 and az_1_2 > slipface_normal:
                direction = 'left'
            elif az_1_0 > slipface_normal and slipface_normal > az_1_2:
                direction = 'right'
            elif slipface_normal > az_1_0 and az_1_0 > az_1_2:
                direction = 'left'
            else:
                print('ERROR: slipface downhill direction azimuth error!')
        
        # operate program for first crestnode
        elif target_crestnode_index == 0:
            # terminology: crestnode 0 is target, crestnode 1 is one ahead
            
            # get azimuth from crestnode 0 to crestnode 1
            az_0_1 = self.calc_azimuth(target_crestnode_index, (target_crestnode_index + 1))
                    
            # get slipface normal from the crestnode (which should be correct)
            slipface_normal = self.crestnodes[target_crestnode_index].slipface_normal
            
            # set conditions based on assumption that slipface normal will be 90 degrees
            # to the end crestnode
            if (slipface_normal - az_0_1) < 90.5 and (slipface_normal - az_0_1) > 89.5:
                direction = 'right'
            elif (az_0_1 - slipface_normal) < 90.5 and (az_0_1 - slipface_normal) > 89.5:
                direction = 'left'
            elif (slipface_normal - az_0_1) < 270.5 and (slipface_normal - az_0_1) > 269.5:
                direction = 'left'
            elif (az_0_1 - slipface_normal) < 270.5 and (az_0_1 - slipface_normal) > 269.5:
                direction = 'right'
            else:
                print('ERROR: slipface downhill direction azimuth error!')
                   
        # operate program for last crestnode
        elif target_crestnode_index == (self.num_crestnodes - 1):
            # terminology: crestnode 1 is target, crestnode 0 is one behind
            
            # get azimuth from crestnode 1 to crestnode 0
            az_1_0 = self.calc_azimuth(target_crestnode_index, (target_crestnode_index - 1))
                    
            # get slipface normal from the crestnode (which should be correct)
            slipface_normal = self.crestnodes[target_crestnode_index].slipface_normal
            
            # set conditions based on assumption that slipface normal will be 90 degrees
            # to the end crestnode
            if (slipface_normal - az_1_0) < 90.5 and (slipface_normal - az_1_0) > 89.5:
                direction = 'left'
            elif (az_1_0 - slipface_normal) < 90.5 and (az_1_0 - slipface_normal) > 89.5:
                direction = 'right'
            elif (slipface_normal - az_1_0) < 270.5 and (slipface_normal - az_1_0) > 269.5:
                direction = 'right'
            elif (az_1_0 - slipface_normal) < 270.5 and (az_1_0 - slipface_normal) > 269.5:
                direction = 'left'
            else:
                print('ERROR: slipface downhill direction azimuth error!')
        
        else:
            print('ERROR: improper target crestnode index call for downhill direction calc')
        
        return direction
  
    def calc_slipface_height(self, crestnode_index, normal_az):
        # method to calculate the intersection distance between a point
        # and the bottom of slipface line feature. Returns 'None' if
        # no intersection is found (indicating we looked in the wrong direction)
        # The algorithm sends out a ray ('line') along the azimuth and direction
        # to look for the intersect. The slipface height is subsequently
        # calculated from the assumption that the slipface holds at the 
        # angle of repose.
        
        # method to calculate the slipface height in initial setup
        # this function calculates slipface height based on the spatial
        # distance between the crest and the bottom of crest feature
        # Note: there is the obvious problem of assuming imagery
        # was orthorectified to the required level of detail. I'm not
        # sure a better method yet without knowing the position of the
        # satellite in space when the image was collected and redoing the
        # orthorectification. For now, I'm going to give this a go.
        
        ray_distance_m = 500.0        # distance in meters to send the ray
        
        # get the origin crestnode position in space
        lon, lat, z = self.crestnodes[crestnode_index].geom.GetPoint()
        
        # set a geopy 'point' object and convert distance
        geopy_point1 = Point(longitude = lon, latitude = lat)
        ray_distance = ray_distance_m / 1000.0
        
        # compute 'vincenty' distance to the new point
        geopy_point2 = VincentyDistance(kilometers = ray_distance).destination(geopy_point1, normal_az)
        
        # construct the ray
        ray = ogr.Geometry(ogr.wkbLineString)
        ray.AddPoint(x = geopy_point1.longitude, y = geopy_point1.latitude, z = z)
        ray.AddPoint(x = geopy_point2.longitude, y = geopy_point2.latitude, z = z)
        
        # evaluate intersection and exit by returning 'None' if no intersect
        if not ray.Intersect(self.bot_geom):
            return None
        else:
            # deal with the case where there is multiple intersections
            ray_inters = ray.Intersection(self.bot_geom)
            num_intersections = ray_inters.GetGeometryCount()
            if num_intersections != 0:
                # ct the closest intersection point by testing them all
                min_test_offsetdist = ray_distance_m
                min_test_index = -1
                for i in range(0, num_intersections):
                    # calculate distance to the intersection
                    test_offset_point = ray_inters.GetGeometryRef(i)
                    test_off_lon, test_off_lat, test_off_z = test_offset_point.GetPoint()
                    geopy_testoffpoint = Point(longitude = test_off_lon, latitude = test_off_lat)
                    test_offsetdist = VincentyDistance(geopy_point1, geopy_testoffpoint).m
                    if test_offsetdist < min_test_offsetdist:
                        min_test_offsetdist = test_offsetdist
                        min_test_index = i
                    
                # ct the closest point based on the min test index variable set above
                offset_point = ray_inters.GetGeometryRef(min_test_index)
                
            else:
                # else set the offset point in the single intersection case
                offset_point = ray.Intersection(self.bot_geom)
            
            # calculate the distance to the intersection, fed from above in both cases
            off_lon, off_lat, off_z = offset_point.GetPoint()
            geopy_offpoint = Point(longitude = off_lon, latitude = off_lat)
            offset_dist = VincentyDistance(geopy_point1, geopy_offpoint).m
            
            # determine the slipface height with trigonometry
            slip_ht = offset_dist * math.tan(math.radians(repose))
            
        return slip_ht
    
    def calc_azimuth(self, from_crestnode_index, to_crestnode_index):
        # convenience function to calculate the azimuth between two crestnodes
        # note: this has to be done with rather elaborate vincenty distance
        #       calculations as latitude and longitude are not equal measures
        
        # Get the locations as lon, lat, z
        lon_A, lat_A, z_A = self.crestnodes[from_crestnode_index].geom.GetPoint()
        lon_B, lat_B, z_B = self.crestnodes[to_crestnode_index].geom.GetPoint()
        
        # Calculate A to B azimuth by setting a fake point to calc x and y
        fakepoint = Point(longitude = lon_A, latitude = lat_B)
        realpoint_A = Point(longitude = lon_A, latitude = lat_A)
        realpoint_B = Point(longitude = lon_B, latitude = lat_B)
        
        distance_x = VincentyDistance(fakepoint, realpoint_B).m
        distance_y = VincentyDistance(fakepoint, realpoint_A).m
        
        # correct distances with positive and negative signs
        if (lon_B - lon_A) < 0.0:
            distance_x = -1.0 * distance_x
        if (lat_B - lat_A) < 0.0:
            distance_y = -1.0 * distance_y
        
        azimuth = math.atan2(distance_x, distance_y)
        azimuth = math.degrees(azimuth)
        
        # correct azimuths to 0-360 scale
        if azimuth < 0.0:
            azimuth = azimuth + 360.0

        return azimuth
    
    def flip_az(self, in_az):
        # utility method to flip an azimuth 180.0 degrees
        out_az = in_az + 180.0
        if out_az > 360.0:
            out_az = out_az - 360.0
        
        return out_az
    
    def advance(self, input_environment, timestep):
        # pass input_input environment and timestep down to each crestnode
        for i in range(0, self.num_crestnodes):
            self.crestnodes[i].advance(input_environment, timestep)

        self.update_geom()           # update the geometry
        return
    
    def calc_curvature(self):
        # method to internally calculate the curvature of the crest
        # and pass the curvature down to each crestnode
        
        # the slipface curvature is the 'gamma' angle from previous
        # work, which denotes a measurement of how curved the slipface
        # point is. This is a measure of the difference between the
        # azimuths of crestlines previous and ahead of the point.
        
        # the first and last node cannot have any slipface curvature
        self.crestnodes[0].slipface_g = 0.0
        self.crestnodes[(self.num_crestnodes - 1)].slipface_g = 0.0
        self.crestnodes[0].b_up = 1.0
        self.crestnodes[(self.num_crestnodes - 1)].b_up = 1.0
        self.crestnodes[0].delivery_gate = 1.0
        self.crestnodes[(self.num_crestnodes - 1)].delivery_gate = 1.0
        
        # update the crestnode objects for the rest of the crestnodes
        for i in range(1, (self.num_crestnodes - 1)):
            # ----------------------------------------------------------------------
            # 1) Calculate slipface gamma based on the angle difference
            #    Positive gamma: concave slipface (e.g., barchan head)
            #    Negative gamma: convex slipface (e.g., parabolic head)
            
            # Calculate the angles from the node point to the previous and next
            angle_P = self.calc_azimuth(i, (i - 1))
            angle_N = self.calc_azimuth(i, (i + 1))
            
            # Calc the difference between the two angles to the left and right
            if angle_P > angle_N:
                l_anglespread = angle_N + 360.0 - angle_P
                r_anglespread = angle_P - angle_N
            else:
                l_anglespread = angle_N - angle_P
                r_anglespread = angle_P + 360.0 - angle_N
            
            if self.downhill_direction == 'left':
                gamma = 180.0 - l_anglespread
            elif self.downhill_direction == 'right':
                gamma = 180.0 - r_anglespread
            else:
                print('ERROR in gamma calculation: need downhill direction')
                
            self.crestnodes[i].slipface_g = gamma   # assign it
            
            # Notes: For the trapezoid area calculation, we set b_up as the distance between
            # the two bisectors on the crestline.
            
            # For the conic area calculation, instead we calculate the 'delivery gate'
            # within each crestnode object. The 'delivery gate' is the cross wind
            # distance between the nearest bisector point, and an adjacent point on
            # the other side of the crest that is equal distance away from the crestnode.
            
            # ----------------------------------------------------------------------
            # 2) Determine the b_up for each crestnode (for trapezoid area calculations)
            
            # Note: we assume plane, just average the latitude and longitude to calculate
            # the bisector points. This is an assumption, but shouldn't make a measurable
            # difference as the scale for this calculation will always be very small.

            # TO DO: fix this portion of the code to use gatepoints as in the conic
            # slipface area calculation method below.
            
            point_P = self.crestnodes[i - 1].geom
            point_target = self.crestnodes[i].geom
            point_N = self.crestnodes[i + 1].geom
            
            # ct the latitude and longitudes from each geometry
            P_lon, P_lat, P_z = point_P.GetPoint()
            target_lon, target_lat, target_z = point_target.GetPoint()
            N_lon, N_lat, N_z = point_N.GetPoint()
            
            # calculate bisectpoint P
            bisectpoint_P_lon = (P_lon + target_lon) / 2.0
            bisectpoint_P_lat = (P_lat + target_lat) / 2.0
            
            # calculate bisectpoint N
            bisectpoint_N_lon = (N_lon + target_lon) / 2.0
            bisectpoint_N_lat = (N_lat + target_lat) / 2.0
            
            # calculate Vincenty distance between the two points to define trapezoid b_up
            geopy_bisect_P = Point(longitude = bisectpoint_P_lon, latitude = bisectpoint_P_lat)
            geopy_bisect_N = Point(longitude = bisectpoint_N_lon, latitude = bisectpoint_N_lat)
            b_up = VincentyDistance(geopy_bisect_P, geopy_bisect_N).m

            # assign to value to the crestnode object
            self.crestnodes[i].b_up = b_up
            
            # ----------------------------------------------------------------------
            # 3) Conic slipface area calculations
            # We first find the closest bisector, this is one gatepoint. The second 
            # gatepoint is then defined as the point which is the same distance away
            # from the node as the first gatepoint, but the opposite direction along
            # the crest. Note this is different than the previous trapezoid calc, which
            # is subject to problems with variable length line segments.
            geopy_target = Point(longitude = target_lon, latitude = target_lat)
            bisect_dist_P = VincentyDistance(geopy_target, geopy_bisect_P).m
            bisect_dist_N = VincentyDistance(geopy_target, geopy_bisect_N).m
            
            # split program flow based on which is closer
            if bisect_dist_P > bisect_dist_N:
                # the 'next' bisect point is closer
                gatepoint_dist = bisect_dist_N
                gatepoint_N_lon = bisectpoint_N_lon
                gatepoint_N_lat = bisectpoint_N_lat
                geopy_gatepoint_N = Point(longitude = gatepoint_N_lon, latitude = gatepoint_N_lat)
                
                # assign the 'previous' gatepoint
                gatepoint_dist_km = gatepoint_dist / 1000.0
                geopy_gatepoint_P = VincentyDistance(kilometers = gatepoint_dist_km).destination(geopy_target, angle_P)
                
                gatepoint_P_lon = geopy_gatepoint_P.longitude
                gatepoint_P_lat = geopy_gatepoint_P.latitude
                                
            else:
                # the 'previous' bisect point is closer
                gatepoint_dist = bisect_dist_P
                gatepoint_P_lon = bisectpoint_P_lon
                gatepoint_P_lat = bisectpoint_P_lat
                geopy_gatepoint_P = Point(longitude = gatepoint_P_lon, latitude = gatepoint_P_lat)
            
                # assign the 'next' gatepoint
                gatepoint_dist_km = gatepoint_dist / 1000.0
                geopy_gatepoint_N = VincentyDistance(kilometers = gatepoint_dist_km).destination(geopy_target, angle_N)
                
                gatepoint_N_lon = geopy_gatepoint_N.longitude
                gatepoint_N_lat = geopy_gatepoint_N.latitude
            
            # now we have the two gatepoints, we can assign the left and right gatepoint
            # variables for each crestnode
            if self.downhill_direction == 'left':
                self.crestnodes[i].gatept_L_lon = gatepoint_P_lon
                self.crestnodes[i].gatept_L_lat = gatepoint_P_lat
                self.crestnodes[i].gatept_R_lon = gatepoint_N_lon
                self.crestnodes[i].gatept_R_lat = gatepoint_N_lat
            elif self.downhill_direction == 'right':
                self.crestnodes[i].gatept_L_lon = gatepoint_N_lon
                self.crestnodes[i].gatept_L_lat = gatepoint_N_lat
                self.crestnodes[i].gatept_R_lon = gatepoint_P_lon
                self.crestnodes[i].gatept_R_lat = gatepoint_P_lat
            else:
                print('ERROR with gatepoint assignment, need downhill direction')
            
            # ----------------------------------------------------------------------
            # Calculate the center of the area circles
            # Previously we had problems using ray_intersects, which seem to be subject
            # to problems with geoid errors and point positioning - we now use a much
            # better system based on geometry.
            
            # We calculate the r_inner based on gamma and the gate distance
            self.crestnodes[i].delivery_gate = VincentyDistance(geopy_gatepoint_P, geopy_gatepoint_N).m

            # First, determine the half_gatedist and absolute value of gamma
            half_gatedist = self.crestnodes[i].delivery_gate / 2.0
            abs_gamma = self.crestnodes[i].slipface_g
            if abs_gamma < 0.0:
                abs_gamma = -1.0 * abs_gamma
            
            # Next, we can calculate the radius to centerpoint
            r_to_center = half_gatedist / math.cos(math.radians((180.0 - abs_gamma) / 2.0))
            
            # We can calculate the positioning of the centerpoint with a point offset
            center_offset = math.sqrt(math.pow(gatepoint_dist, 2.0) + math.pow(r_to_center, 2.0))
            center_offset_km = center_offset / 1000.0
            if self.crestnodes[i].slipface_g > 0.0:
                center_offset_az = self.crestnodes[i].slipface_normal
            else:
                center_offset_az = self.flip_az(self.crestnodes[i].slipface_normal)
            geopy_center_pt = VincentyDistance(kilometers = center_offset_km).destination(geopy_target, center_offset_az)
            
            # and . . assign it to the crestnode object for posterity
            self.crestnodes[i].curve_centerpoint_lon = geopy_center_pt.longitude
            self.crestnodes[i].curve_centerpoint_lat = geopy_center_pt.latitude
            
            # We can also calculate the open angle, as it is equal to abs_gamma
            self.crestnodes[i].open_angle = abs_gamma
            
            # A quick check - the trig breaks down if we go too far out and precision errors
            # creep in to play is we have centerpoints that are too far away.
            max_center_offset_dist = 5000.0
            if center_offset > max_center_offset_dist:
                
                # here, the slipface is basically flat.
                # assign the centerpoint to be None, to signal to the crestnode advance function
                # to calculate the slipface area as a plane.
                self.crestnodes[i].curve_centerpoint_lon = None
                self.crestnodes[i].curve_centerpoint_lat = None
                self.crestnodes[i].r_inner = None
                self.crestnodes[i].r_outer = None
                self.crestnodes[i].open_angle = None

            else:
                # Else, proceed normally and assign the inner and outer radii . .            
                # to calculate the inner and outer radii, we need to assess which way is down
                # if the slipface gamma is positive, this means that the intersect point is downhill
                # and the outer radius is at the top of the slipface, the inner radius is at the bottom
                if self.crestnodes[i].slipface_g > 0.0:
                    # concave slipface
                    self.crestnodes[i].r_outer = r_to_center
                    r_inner = self.crestnodes[i].r_outer - self.crestnodes[i].w
                    
                    # check to see if we have an intersection before the bottom
                    if r_inner < 0.0:
                        r_inner = 0.0       # correct it to zero
                        
                    self.crestnodes[i].r_inner = r_inner
                    
                elif self.crestnodes[i].slipface_g < 0.0:
                    # convex slipface
                    self.crestnodes[i].r_inner = r_to_center
                    self.crestnodes[i].r_outer = self.crestnodes[i].r_inner + self.crestnodes[i].w
       
        return
               
    def update_geom(self):
        # wrapper method to update all node geometric attributes
        # and internally update the list of point geometries to
        # a linestring feature.
        
        # update internal geom object
        for i in range(0, self.num_crestnodes):
            pt_lon, pt_lat, pt_z = self.crestnodes[i].geom.GetPoint()
            self.geom.SetPoint(i, x = pt_lon, y = pt_lat, z = pt_z)
        
        # update slipface normal, and calculate heights
        self.calc_normals_heights()
        
        # set the downhill direction
        if self.downhill_direction == None:
            self.set_downhill_direction()
        
        # once we have the downhill direction set, calculate the curvature
        self.calc_curvature()

        return
    
    def debug(self, input_environment, timestep):
        # method to push debugging information for interrogation
        # creates a kml file with the crestnodes and relevant information
        # each crest is it's own kml file with a series of node points
        
        # first thing to do is to run a fake advance calculation over
        # the crestnodes to update all the internals (e.g., flux, etc.)
        # These internals are updated during the advance call, and represent
        # the values that resulted in the present positioning of the crestnodes
        # By running a fake call to advance that calculates everything,
        # except doesn't move the nodes, we can produce some values to
        # interrogate the calculations.
        
        for i in range(0, self.num_crestnodes):
            self.crestnodes[i].bluff_advance(input_environment, timestep)
        
        if extra_detailed_debugging:
            print('WARNING: outputting extra detailed debugging data')
        
        # check to see if the debug directory exists
        if not os.path.isdir('debug'):
            os.mkdir('debug')
            
        # enter the directory to place the file
        base_dir = os.getcwd()
        os.chdir(os.path.join(os.getcwd(), 'debug'))
        
        # construct the filename: debug_time_crest_crestid.kml
        time_string = str(input_environment.get_time(timestep))
        debug_filename = 'debug_' + time_string + '_crest_' + str(self.id) + '.kml'
        
        out_drv = ogr.GetDriverByName('KML')
        out_data = out_drv.CreateDataSource(debug_filename)
        out_layer = out_data.CreateLayer('debug', geom_type=ogr.wkbPoint)
        
        # create the layers
        out_layer.CreateField(ogr.FieldDefn('timestep', ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn('crest_id', ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn('crestnode_id', ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn('slipface_normal', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('slipface_g', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('veg', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('q_vol', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('q_az', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('move_az', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('w', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('phi', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('b_up', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('b_down', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('dh_dt_slipface', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('veg_threshold', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('movedist', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('flux_down_slipface', ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn('sed_delivered', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('slipface_area', ogr.OFTReal))

        # conic area variables
        out_layer.CreateField(ogr.FieldDefn('gatept_L_lon', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('gatept_L_lat', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('gatept_R_lon', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('gatept_R_lat', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('r_inner', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('r_outer', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('curve_centerpoint_lon', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('curve_centerpoint_lat', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('delivery_gate', ogr.OFTReal))
        out_layer.CreateField(ogr.FieldDefn('open_angle', ogr.OFTReal))

        out_featuredef = out_layer.GetLayerDefn()
        
        # loop through the crestnodes and store the results
        for i in range(0, self.num_crestnodes):
            # create feature and define the fields
            out_feat = ogr.Feature(out_featuredef)
            out_feat.SetGeometry(self.crestnodes[i].geom)

            # set the fields
            out_feat.SetField('timestep', timestep)
            out_feat.SetField('crest_id', self.id)
            out_feat.SetField('crestnode_id', i)
            out_feat.SetField('slipface_normal', self.crestnodes[i].slipface_normal)
            out_feat.SetField('slipface_g', self.crestnodes[i].slipface_g)
            out_feat.SetField('veg', self.crestnodes[i].veg)
            out_feat.SetField('b_up', self.crestnodes[i].b_up)
            out_feat.SetField('q_vol', self.crestnodes[i].q_vol)
            out_feat.SetField('q_az', self.crestnodes[i].q_az)
            out_feat.SetField('veg_threshold', self.crestnodes[i].veg_threshold)
            out_feat.SetField('w', self.crestnodes[i].w)
            out_feat.SetField('phi', self.crestnodes[i].phi)
            out_feat.SetField('b_down', self.crestnodes[i].b_down)
            out_feat.SetField('dh_dt_slipface', self.crestnodes[i].dh_dt_slipface)
            out_feat.SetField('movedist', self.crestnodes[i].movedist)
            out_feat.SetField('flux_down_slipface', str(self.crestnodes[i].flux_down_slipface))
            out_feat.SetField('sed_delivered', self.crestnodes[i].sed_delivered)
            out_feat.SetField('slipface_area', self.crestnodes[i].slipface_area)
            out_feat.SetField('move_az', self.crestnodes[i].move_az)
            
            # conic area calculation variables
            out_feat.SetField('gatept_L_lon', self.crestnodes[i].gatept_L_lon)
            out_feat.SetField('gatept_L_lat', self.crestnodes[i].gatept_L_lat)
            out_feat.SetField('gatept_R_lon', self.crestnodes[i].gatept_R_lon)
            out_feat.SetField('gatept_R_lat', self.crestnodes[i].gatept_R_lat)
            out_feat.SetField('r_inner', self.crestnodes[i].r_inner)
            out_feat.SetField('r_outer', self.crestnodes[i].r_outer)
            out_feat.SetField('curve_centerpoint_lon', self.crestnodes[i].curve_centerpoint_lon)
            out_feat.SetField('curve_centerpoint_lat', self.crestnodes[i].curve_centerpoint_lat)
            out_feat.SetField('delivery_gate', self.crestnodes[i].delivery_gate)
            out_feat.SetField('open_angle', self.crestnodes[i].open_angle)
            
            out_layer.CreateFeature(out_feat)
            
            #  detailed debugging information (only run with one crest, and one timestep!)
            if extra_detailed_debugging and i > 0 and i < (self.num_crestnodes - 1):
                
                # push center of curve points
                self.simple_kmlpointpush(filename = 'center_' + str(i) + '.kml',
                                         lon = self.crestnodes[i].curve_centerpoint_lon,
                                         lat = self.crestnodes[i].curve_centerpoint_lat
                                         )
                # push L side gatepoint
                self.simple_kmlpointpush(filename = 'gatepoint_L_' + str(i) + '.kml',
                                         lon = self.crestnodes[i].gatept_L_lon,
                                         lat = self.crestnodes[i].gatept_L_lat
                                         )

                # push R side gatepoint
                self.simple_kmlpointpush(filename = 'gatepoint_R_' + str(i) + '.kml',
                                         lon = self.crestnodes[i].gatept_R_lon,
                                         lat = self.crestnodes[i].gatept_R_lat
                                         )
                
        out_data.Destroy()
        
        # back to original directory        
        os.chdir(base_dir)
        
        return

    def simple_kmlpointpush(self, filename, lon, lat):
        # utility method: simple method to push out a kml file of a point for debugging
        
        drv = ogr.GetDriverByName('KML')
        data = drv.CreateDataSource(filename)
        layer = data.CreateLayer('debug', geom_type=ogr.wkbPoint)
                        
        featuredef = layer.GetLayerDefn()              
        feat = ogr.Feature(featuredef)
        geom = ogr.Geometry(ogr.wkbPoint)
        geom.AddPoint(x = lon, y = lat, z = 0.0)
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        data.Destroy()
        
        return




        