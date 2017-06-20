# This file is part of BrownianDynamics. 
# Copyright (C) 2017 Manuel Luitz <manuel.luitz@gmail.com>
# Copyright (C) 2017 Fabian Zeller <fabian.zeller11@gmail.com>
#
# BrownianDynamics is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# BrownianDynamics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with BrownianDynamics. If not, see <http://www.gnu.org/licenses/>.
# 
import random
import copy
from numpy.linalg import norm

import BD_propagate as p
import BD_pqr as pqr

def trajectory(args):
    i, ligand_prototype, receptor, grid, sig_relative_trans, sig_ligand_rot, sig_receptor_rot, CONFIG, queue = args
    # Preparation -------------------------------------------------------------
    
    # initiate a separate instance of Random for every trajectory,
    # so that in the case of threading every thread gets different random numbers
    ran = random.Random()
    ran.jumpahead(CONFIG.NUMBER_TRAJECTORIES)
    
    # create an own copy of the propagated object for every trajectory run
    ligand = copy.deepcopy(ligand_prototype)
    
    # move it to a random position on a sphere of radius STARTING_RADIUS
    # with random orientation and update grid indices
    p.starting_position(ligand, CONFIG, ran)
    ligand.updateGridIndices(grid)
    
    # check initial state of the system
    if p.collision(ligand, receptor, grid):
        print 'error: starting position collides'
    state_tmp = p.state(ligand, receptor, CONFIG)
    if not state_tmp ==2:
        'error: starting position not between rmax and receptor'

    # Propagation Loop --------------------------------------------------------
    while state_tmp == 2:
        ligand.savePreviousPosition()
        
        # Adjust step size for if ligand and receptor are separated
        if norm(ligand.center - receptor.center) > ligand.rmax + receptor.rmax:
            sig_factor = 3.33
        else:
            sig_factor = 1.
        
        # Propagation 
        p.propagate_trans(ligand, sig_relative_trans*sig_factor, ran)
        p.propagate_rot_receptor(ligand, sig_receptor_rot*sig_factor, ran)
        p.propagate_rot_ligand(ligand, sig_ligand_rot*sig_factor, ran)
        ligand.updateGridIndices(grid)
        
        
        if p.collision(ligand, receptor, grid):
            ligand.resetToPreviousPosition()
            continue
        state_tmp = p.state(ligand, receptor, CONFIG)
        
        if state_tmp == 1:
            pqr.write_molecules_to_pqr(ligand, 'target_'+str(i), CONFIG.LIGAND_PQR)

    # -------------------------------------------------------------------------
    del args
    if not queue == 0:
        queue.put(state_tmp)
    return state_tmp
