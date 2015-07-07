# Crestline simulation tool
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

import os
from dunefield import *         # import dunefield class
from environment import *       # import environment class
from xl_integration import *    # excel integration class

from joblib import Parallel, delayed

def welcome():
    print('------------------------------------------------------------------------')
    print('Welcome to DAPPER: the Dune Advance Program for Predicting Earth Reactivations')
    print('Tom Barchyn, Chris Hugenholtz, contact: tbarchyn@gmail.com')
    print('University of Calgary, Calgary, AB, Canada')
    print('This program has no warranty! The authors are not responsible for outputs')
    print('Note: this program only works on Earth (update datum for other planets)')
    print('------------------------------------------------------------------------')
    return

def run():
    # Run function, to be called below
    # set up the excel link
    # Sort out the ridiculous problems with working directory
    # that are created by xlwings
    exec_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(exec_dir)    
    
    # set up the links
    xl_filename = os.path.join(os.getcwd(), 'dapper_control.xlsm')
    xl_control = xl_link(xl_filename)
    
    # grab parameters from the xl file
    input_topkml = xl_control.topkml
    input_bottomkml = xl_control.bottomkml
    outputkml = xl_control.oputkml
    yearlist = xl_control.grab_time()
    oputschedule = xl_control.grab_oputschedule()
    debug_dapper = xl_control.debug
    
    welcome()

    # (1) initialize environment for simulation
    t_env = environment(xl_control)
    
    # (2) initialize dune field from input
    t_dunefield = dunefield(input_topkml, input_bottomkml)
    
    # (3) run the dune field forward in time
    for t in range(0, t_env.timesteps):
        t_dunefield.advance(t_env, t)
        
        # call debug code to push debug outputs
        if debug_dapper:
            t_dunefield.debug(t_env, t)
        
        # check output schedule to see if its time to output interim file
        if oputschedule[t]:
            interim_oputfile = xl_control.interim_prefix + str(yearlist[t]) + '.kml'
            interim_oputfile = interim_oputfile
            t_dunefield.output(interim_oputfile)
        
        # message the user
        user_message = 'Timestep completed: ' + str(yearlist[t])
        print(user_message)
        xl_control.message(user_message)
    
    # (4) output final dune field crest position 
    t_dunefield.output(outputkml)
    print ('Simulation finished!')
    xl_control.message('Simulation finished . . ')
    
    return
    
if __name__ == '__main__':
    run()
    











