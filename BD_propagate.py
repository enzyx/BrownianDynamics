import numpy
import math
import BD_geometry as g
import BD_diffusion as d

from numpy.linalg import norm

class Propagator(object):
    def __init__(self, ligand, receptor, CONFIG):
        self.sigma_trans_relative  = d.sigma_relative_trans(ligand.rgyr, receptor.rgyr, CONFIG)  
        self.sigma_rot_ligand      = d.sigma_rot(ligand.rgyr, CONFIG) 
        self.sigma_rot_receptor    = d.sigma_rot(receptor.rgyr, CONFIG)
        self.sigma_trans_relative2 = self.sigma_trans_relative * math.sqrt(CONFIG.DT_BEYOND_INTERA_RADIUS / CONFIG.DT)
        self.sigma_rot_ligand2     = self.sigma_rot_ligand     * math.sqrt(CONFIG.DT_BEYOND_INTERA_RADIUS / CONFIG.DT)
        self.sigma_rot_receptor2   = self.sigma_rot_receptor   * math.sqrt(CONFIG.DT_BEYOND_INTERA_RADIUS / CONFIG.DT)
        
    def propagate_trans(self, ligand, ran, interacting, F):
        """
        calculate the RELATIVE translational diffusion 
        of receptor and ligand,
        apply it to ligand atom positions
        """
        if interacting == True:
            dR = numpy.array([ran.gauss(0., self.sigma_trans_relative),
                              ran.gauss(0., self.sigma_trans_relative),
                              ran.gauss(0., self.sigma_trans_relative)])
            dR += F * self.sigma_trans_relative
        else: 
            dR = numpy.array([ran.gauss(0., self.sigma_trans_relative2),
                              ran.gauss(0., self.sigma_trans_relative2),
                              ran.gauss(0., self.sigma_trans_relative2)]) 
                      
        g.translate(ligand, dR)
        
        
    def propagate_rot_receptor(self, ligand, ran, interacting, T_rec):
        """
        calculate the rotational diffusion of the
        receptor and apply it to ligand atom positions
        """ 
        if interacting == True:     
            dO = numpy.array([ran.gauss(0., self.sigma_rot_receptor),
                              ran.gauss(0., self.sigma_rot_receptor),
                              ran.gauss(0., self.sigma_rot_receptor)])
            dO += T_rec * self.sigma_rot_receptor
        else:
            dO = numpy.array([ran.gauss(0., self.sigma_rot_receptor2),
                              ran.gauss(0., self.sigma_rot_receptor2),
                              ran.gauss(0., self.sigma_rot_receptor2)])   
                 
        g.rotate(ligand, dO)
    
    def propagate_rot_ligand(self, ligand, ran, interacting, T_lig):
        """
        calculate the rotational diffusion of the
        ligand and apply it to ligand atom positions
        """      
        if interacting == True: 
            dO = numpy.array([ran.gauss(0., self.sigma_rot_ligand),
                              ran.gauss(0., self.sigma_rot_ligand),
                              ran.gauss(0., self.sigma_rot_ligand)])
            dO += T_lig * self.sigma_rot_ligand
        else:
            dO = numpy.array([ran.gauss(0., self.sigma_rot_ligand2),
                              ran.gauss(0., self.sigma_rot_ligand2),
                              ran.gauss(0., self.sigma_rot_ligand2)])
            
        # go to ligand origin, rotate, go back to original position
        g.rotate(ligand, dO, True)
    
    
    def state(self, ligand, receptor, CONFIG):
        
        ligand_radius           = norm(ligand.center)
        if ligand_radius > CONFIG.MAXIMUM_RADIUS:
            return 0
        else:
            reaction_coordinate = self.getReactionCoordinate(ligand, receptor, CONFIG)
            #print reaction_coordinate
            if reaction_coordinate < CONFIG.COORD_THRESHOLD:
                return 1
            else:
                return 2
    
    def getReactionCoordinate(self, ligand, receptor, CONFIG):
        # atom pair distance rmsd
        if CONFIG.REACTION_COORD_TYPE == 'drmsd':
            drmsd   = 0.0
            N_pairs = len(CONFIG.LIGAND_COORD_ATOMS)
            for i in range(N_pairs):
                drmsd += norm(ligand.R[CONFIG.LIGAND_COORD_ATOMS[i]] - \
                              receptor.R[CONFIG.RECEPTOR_COORD_ATOMS[i]])**2
            return (drmsd/N_pairs)**0.5
        # com is default
        else:
            ligand_coord        = g.center(ligand.R, CONFIG.LIGAND_COORD_ATOMS) 
            receptor_coord      = g.center(receptor.R, CONFIG.RECEPTOR_COORD_ATOMS)
            return norm(ligand_coord - receptor_coord)
        

            
    def starting_position(self, ligand, CONFIG, ran):
        
        def randomSinusDistribution():
            y = 1.
            x = 0.
            while y > math.sin(x):
                x = ran.uniform(0, 2*math.pi)
                y = ran.uniform(0., 1.)
            return x
        
        g.translate(ligand, - ligand.center) 
        
        # Move ligand to starting radius on z-axis
        dR = numpy.array([0., 0., CONFIG.STARTING_RADIUS])
        g.translate(ligand, dR)
         
        # random point on sphere with radius starting_radius
        dO = numpy.matrix([[0.],
                           [0.],
                           [ran.uniform(0, 2*math.pi)]])
        g.rotate(ligand, dO)
        dO = numpy.matrix([[0.],
                           [randomSinusDistribution()],
                           [0.]])
        g.rotate(ligand, dO)
        
        dO = numpy.matrix([[ran.uniform(0, 2*math.pi)],
                          [0.],
                          [0.]])
        g.rotate(ligand, dO)
         
        
       
       
       
#         dR = numpy.array([31., 0., 0.])
#         g.translate(ligand, dR)
#  
#         dO = numpy.matrix([[-0.7*math.pi],
#                            [0.2*math.pi],
#                             [0.13*math.pi]])
#         g.rotate(ligand, dO)
#          
#         dO = numpy.matrix([[ran.uniform(-0.03*math.pi, 0.03*math.pi)],
#                             [ran.uniform(-0.03*math.pi, 0.03*math.pi)],
#                             [ran.uniform(-0.03*math.pi, 0.03*math.pi)]])
#         g.rotate(ligand, dO)
         
        
        
        # random orientation of ligand
        
        dO = numpy.matrix([[0.],
                            [0.],
                            [ran.uniform(0, 2*math.pi)]])
        g.rotate(ligand, dO, True)
        dO = numpy.matrix([[0.],
                            [randomSinusDistribution()],
                            [0.]])
        g.rotate(ligand, dO, True)
       
        dO = numpy.matrix([[ran.uniform(0, 2*math.pi)],
                           [0.],
                           [0.]])
        g.rotate(ligand, dO, True)  
        
        ligand.updateCenter()
    
        
