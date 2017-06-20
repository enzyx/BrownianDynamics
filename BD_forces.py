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
import numpy
from numpy.linalg import norm


def calcForceAndTorque(ligand, receptor, grid, CONFIG):
    """
    look up electrostatic field in the grid
    calculate force and torques
    """
    F       = numpy.zeros([3], float)
    T_lig   = numpy.zeros([3], float)
    T_rec   = numpy.zeros([3], float)
    
    for i in range(ligand.N_atoms):
        # Is the atom within the grid?
        if ligand.grid_indices[i,0] > 0 and ligand.grid_indices[i,1] > 0 and ligand.grid_indices[i,2] > 0:
            # FORCE
            dF = - grid.grid[ligand.grid_indices[i,0], 
                             ligand.grid_indices[i,1], 
                             ligand.grid_indices[i,2],
                             :3] * (ligand.q[i])
            F += dF
            
            # TORQUE LIGAND
            r = ligand.R[i] - ligand.center
            T_lig[0] += r[1]*dF[2] - r[2]*dF[1]
            T_lig[1] += r[2]*dF[0] - r[0]*dF[2]
            T_lig[2] += r[0]*dF[1] - r[1]*dF[0]
            
            # TORQUE RECEPTOR
            T_rec[0] -= ligand.R[i,1]*dF[2] - ligand.R[i,2]*dF[1]
            T_rec[1] -= ligand.R[i,2]*dF[0] - ligand.R[i,0]*dF[2]
            T_rec[2] -= ligand.R[i,0]*dF[1] - ligand.R[i,1]*dF[0]
                        
    return F, T_lig, T_rec

def collision(ligand, receptor, grid):
    # Can the molecules collide?
    if norm(ligand.center) - ligand.rmax - receptor.rmax < 0.:
    # Yes: go through all ligand atoms
        for i in range(ligand.N_atoms):
            # Is the atom within the grid?
            if ligand.grid_indices[i,0] > 0 and ligand.grid_indices[i,1] > 0 and ligand.grid_indices[i,2] > 0:
                # Does the atom collide?
                if grid.grid[ligand.grid_indices[i,0], 
                             ligand.grid_indices[i,1], 
                             ligand.grid_indices[i,2],
                             3] == True:
                    return  True
    return False       

def interacting(ligand, CONFIG):
    if norm(ligand.center) >= CONFIG.INTERACTION_RADIUS:
        return False
    else:
        return True
