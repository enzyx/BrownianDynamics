import numpy
import math
import copy
import BD_geometry as g
import BD_diffusion as d
import BD_forces as f
import BD_constants as c
import BD_pqr as pqr
from math import cos, sin
import sys
from numpy import unravel_index

class Ligand(object):
    """
    Ligand Class
    """
    def __init__(self, ligand_pqr, grid):
        
        self.R, self.q          = pqr.load_molecule_from_pqr(ligand_pqr) # N_atoms*3 float 
        self.N_atoms            = len(self.R)                                   # integer
        self.center             = g.center(self.R)                              # 1*3 array float
        self.grid_indices       = numpy.zeros([self.N_atoms,3], int)            # 1*3 array integer
        
        # Initialize previous position properties
        self.savePreviousPosition()

        # Calculate static geometric properties        
        self.rgyr    = g.rgyr(self)
        self.rmax    = g.rmax(self)
        
    def savePreviousPosition(self):
        self.R_prev             = numpy.copy(self.R)
        self.center_prev        = numpy.copy(self.center) 
        self.grid_indices_prev  = numpy.copy(self.grid_indices)

    def resetToPreviousPosition(self):
        self.center             = self.center_prev
        self.R                  = self.R_prev
        self.grid_indices       = self.grid_indices_prev
        
    def updateGridIndices(self, grid):
        for i in range(self.N_atoms):
            for d in range(3):
                index = int( (self.R[i,d] - grid.R_min[d]) / grid.grid_spacing)
                if index < 0.:
                    index = -1
                elif index > grid.shape[d] -1:
                    index = -2
                self.grid_indices[i,d] =  index
    
    def updateCenter(self):
        self.center = g.center(self.R)
                

class Receptor(object):
    """
    Receptor Class
    """
    def __init__(self, CONFIG):

        self.R, self.q          = pqr.load_molecule_from_pqr(CONFIG.RECEPTOR_PQR)
        self.N_atoms            = len(self.R)
        self.center             = g.center(self.R)
        self.center_original    = copy.deepcopy(self.center)
        
        # Center molecule 
        g.translate(self, -self.center)

        # Calculate static geometric properties 
        self.rgyr    = g.rgyr(self)
        self.rmax    = g.rmax(self)

class Grid(object):
    """
    Grid class that contains the receptor atoms divided into a grid
    """
    def __init__(self, receptor, CONFIG):
        # LOAD POTENTIAL, COMPUTE GRADIENT
        potential_flattened = numpy.loadtxt(CONFIG.POTENTIAL_FILE).flatten()
        ps = numpy.array(CONFIG.POTENTIAL_GRID_SHAPE)
        pd = numpy.array(CONFIG.POTENTIAL_GRID_DELTA)
        po = numpy.array(CONFIG.POTENTIAL_GRID_ORIGIN) - numpy.array(receptor.center_original)
        potential = numpy.zeros( CONFIG.POTENTIAL_GRID_SHAPE, float)
        
               
        index = numpy.zeros([3], int)
        for i in range(len(potential_flattened)):
            index[0]  = int(i / ps[1]/ps[2])
            index[1]  = int((i - index[0] * ps[1]*ps[2]) / ps[2])
            index[2]  = int((i - index[0] * ps[1]*ps[2] - index[1] * ps[2]))
            potential[index[0], index[1], index[2]] = potential_flattened[i]
        
        
        gs = copy.deepcopy(ps)
        gs = numpy.append(gs, 3)
        
        gradient = numpy.zeros(gs, float)
        for x in range(gs[0]-1):
            for y in range(gs[1]-1):
                for z in range(gs[2]-1):
                    gradient[x,y,z,0] = (potential[x+1,y,z] - potential[x,y,z]) / pd[0]
                    gradient[x,y,z,1] = (potential[x,y+1,z] - potential[x,y,z]) / pd[1]
                    gradient[x,y,z,2] = (potential[x,y,z+1] - potential[x,y,z]) / pd[2]
                    
       
        # SET UP GRID
        self.grid_spacing   = CONFIG.GRID_SPACING
        self.imax_surr      = int(1. * CONFIG.COLLISION_RADIUS / CONFIG.GRID_SPACING)
        self.imax_surr2     = self.imax_surr**2
        self.R_min          = po
        self.R_max          = po + pd*ps

        self.shape          = numpy.array( map(int, (self.R_max - self.R_min) / self.grid_spacing) ) + 1
        sys.stdout.write('   grid shape: {}\n'.format(self.shape))
        sys.stdout.flush()
        self.grid           = numpy.zeros(numpy.append(self.shape, 4), float)
        
        
        # FILL IN FORCE
        index = numpy.zeros([3], int)
        sys.stdout.write('   Filling in forces...\n')
        sys.stdout.flush()
        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                for z in range(self.shape[2]):
                    index[0] = round(x* self.grid_spacing / pd[0])
                    index[1] = round(y* self.grid_spacing / pd[1]) 
                    index[2] = round(z* self.grid_spacing / pd[2])
                    for d in range(3):
                        if index[d]== gs[d]:
                            index[d] -= 1
                    self.grid[x,y,z, 0] = gradient[index[0], index[1], index[2],0]
                    self.grid[x,y,z, 1] = gradient[index[0], index[1], index[2],1]
                    self.grid[x,y,z, 2] = gradient[index[0], index[1], index[2],2]
        

       
        # FILL IN RECEPTOR ATOMS
        index = numpy.zeros(3, int)
        for i in range(len(receptor.R)):
            sys.stdout.write('\r   inserting atom {:6d}/{:6d}'.format(i+1, len(receptor.R)))
            sys.stdout.flush()
            for d in range(3):
                index[d] = int( (receptor.R[i,d] - self.R_min[d]) / self.grid_spacing)
            for dx in range(-self.imax_surr, self.imax_surr):
                for dy in range(-self.imax_surr, self.imax_surr):
                    for dz in range(-self.imax_surr, self.imax_surr):
                        if dx**2+dy**2+dz**2 < self.imax_surr2:
                            if 0.<index[0]+dx<self.shape[0] and 0.<index[1]+dy<self.shape[1] and 0.<index[2]+dz<self.shape[2]:
                                self.grid[index[0]+dx, index[1]+dy, index[2]+dz, 3] = 1
                
        #numpy.savetxt('receptor', self.grid[:,300,:,3])
        #numpy.savetxt('gradx', self.grid[:,300,:,0]/100)
        #numpy.savetxt('grady', self.grid[:,300,:,1]/100)
        #numpy.savetxt('gradz', self.grid[:,300,:,2]/100)
        
        
        sys.stdout.write('\n')
        sys.stdout.flush()  