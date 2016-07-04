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
