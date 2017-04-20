import ConfigParser
from numpy.linalg import norm
import sys,os,glob
import datetime
import BD_constants as c
import math
from BD_colors import Colors as colors

class Configuration(object):
    """
    contains the configuration parameterse
    """
    def __init__(self, CONFIGFILE):
        if not os.path.isfile(CONFIGFILE):
            print "{}Error:{} Could not open config file: {}".format(colors.FAIL, colors.ENDC, CONFIGFILE)
            sys.exit(-1)

        config = ConfigParser.ConfigParser()
        config.read(CONFIGFILE)
        
        # SYSTEM PROPERTIES
        self.STARTING_RADIUS         = float(       config.get('system', 'starting-radius') )
        self.MAXIMUM_RADIUS          = float(       config.get('system', 'maximum-radius') )
        self.T                       = float(       config.get('system', 'temperature') ) 
        self.COLLISION_RADIUS        = float(       config.get('system', 'collision-radius') )
        self.LIGAND_PQRS             =              self.getLigandPqrs(config)
        self.RECEPTOR_PQR            =              config.get('system', 'receptor.pqr-file')
        self.POTENTIAL_FILE          =              config.get('system', 'potential.dx-file')
        
        # GRID 
        self.GRID_SPACING            = float(       config.get('grid', 'grid-spacing') )
        
        # PICKLE
        self.LOAD_MOLS_FROM_PICKLE   =              config.get('pickle', 'load') 
        self.SAVE_MOLS_TO_PICKLE     =              config.get('pickle', 'save')
        
        # REACTION COORDINATE
        self.LIGAND_COORD_ATOMS      = map(int,     config.get('reaction-distance', 'ligand-atoms').split() )
        self.RECEPTOR_COORD_ATOMS    = map(int,     config.get('reaction-distance', 'receptor-atoms').split() )
        self.COORD_THRESHOLD         = float(       config.get('reaction-distance', 'threshold') )
        self.REACTION_COORD_TYPE     =              config.get('reaction-distance', 'coordinate-type').lower() # possible are: com, drmsd
        self.REFERENCE_DISTANCES     = map(float,   config.get('reaction-distance', 'reference-distances').split() )
        
        # SIMULATION PARAMETERS
        self.RANDOM_START_POSITIONS  = bool(        config.getboolean('simulation', 'random-start-positions') )
        self.DT                      = float(       config.get('simulation', 'integration-time-step') )
        self.DT_BEYOND_INTERA_RADIUS = float(       config.get('simulation', 'integration-time-step-beyond-ia-radius') )
        self.INCLUDE_ELECTROSTATICS  = bool(        config.getboolean('simulation', 'include-electrostatics') )
        self.INTERACTION_RADIUS      = float(       config.get('simulation', 'interaction-radius'))
        self.NUMBER_TRAJECTORIES     = int(         config.get('simulation', 'number-of-trajectories') )
        self.NUMBER_THREADS          = int(         config.get('simulation', 'number-of-threads') )
        self.OUTPUT_DIRECTORY        =              config.get('simulation', 'output-directory')
        self.SAVE_START_COORDS       = bool(        config.getboolean('simulation', 'save-start-coords') )

    def getLigandPqrs(self, config):
        ligand_pqrs = []
        for ligand_pattern in config.get('system', 'ligand.pqr-files').split():
            for ligand_pqr in glob.glob(ligand_pattern):
                ligand_pqrs.append(ligand_pqr)
        return ligand_pqrs
                
def checkInputConsistency(ligand_prototype, receptor, grid, propagator, CONFIG):
    
    if ligand_prototype.rmax + receptor.rmax > CONFIG.INTERACTION_RADIUS:
        print '{ws}WARNING:{we} ligand radius {rl} + receptor radius {rr} > interaction radius {ri}'.format(
                    rl=ligand_prototype.rmax, rr=receptor.rmax, ri=CONFIG.INTERACTION_RADIUS, ws=colors.WARNING, we=colors.ENDC)
        sys.exit()
    if CONFIG.COLLISION_RADIUS < CONFIG.GRID_SPACING:
        print '{ws}WARNING:{we} collision radius < grid resolution'.format(ws=colors.WARNING, we=colors.ENDC)
        sys.exit()
    if CONFIG.REACTION_COORD_TYPE == 'drmsd' and len(CONFIG.LIGAND_COORD_ATOMS) != len(CONFIG.RECEPTOR_COORD_ATOMS) != len(CONFIG.REFERENCE_DISTANCES):
        print '{ws}WARNING:{we} number of ligand and receptor atoms and reference distances'.format(ws=colors.WARNING, we=colors.ENDC)
        print '         need to be equal for drmsd interaction type'
        sys.exit()
        
def printSimulationInfo(ligand_prototype, receptor, grid, propagator, CONFIG):
    
    M_ligand = 22*13*1.6e-27
    
    for output in [sys.stdout, open(CONFIG.OUTPUT_DIRECTORY+'/BD.info','w')]:
        output.write(' Info:\n')
        output.write('  - start time:                {}\n'.format(datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
        output.write('  - Ligand:\n') 
        output.write('      number of atoms:         {:d}\n'.format(ligand_prototype.N_atoms))
        output.write('      total charge:            {:3.2f}\n'.format(sum(ligand_prototype.q)))
        output.write('      radius of gyration:      {:3.2f}\n'.format(ligand_prototype.rgyr))
        output.write('      sigma rotation:          {:3.2f}\n'.format(propagator.sigma_rot_ligand))
        output.write('      sigma translation (rel.) {:3.2f}\n'.format(propagator.sigma_trans_relative))
        output.write('  - Receptor:\n') 
        output.write('      number of atoms:         {:d}\n'.format(receptor.N_atoms))
        output.write('      radius of gyration:      {:3.2f}\n'.format(receptor.rgyr))
        output.write('      sigma rotation    :      {:3.2f}\n'.format(propagator.sigma_rot_receptor))
        output.write('  - Interaction:\n')
        output.write('      reaction coordinate:     {}\n'.format(CONFIG.REACTION_COORD_TYPE))
        output.write('      atoms ligand:            {}\n'.format(CONFIG.LIGAND_COORD_ATOMS))
        output.write('      atoms receptor:          {}\n'.format(CONFIG.RECEPTOR_COORD_ATOMS))
        output.write('      reaction threshold:      {:3.2f}\n'.format(CONFIG.COORD_THRESHOLD))
        output.write('  - Grid:\n')    
        output.write('      grid spacing:            {:3.2f}\n'.format(grid.grid_spacing))    
        output.write('      grid shape:              {}\n'.format(grid.shape))
        output.write('  - Simulation:\n')  
        output.write('      number trajectories:     {:d}\n'.format(CONFIG.NUMBER_TRAJECTORIES))
        output.write('      starting radius:         {:3.2f}\n'.format(CONFIG.STARTING_RADIUS))
        output.write('      max. radius:             {:3.2f}\n'.format(CONFIG.MAXIMUM_RADIUS))  
        output.write('      number start positions:  {}\n'.format(len(CONFIG.LIGAND_PQRS)))
        output.write('      random start positions:  {}\n'.format(CONFIG.RANDOM_START_POSITIONS))  
        output.write('      collision radius:        {:3.2f}\n'.format(CONFIG.COLLISION_RADIUS))    
        output.write('      electrostatics:          {}\n'.format(CONFIG.INCLUDE_ELECTROSTATICS))
        output.write('      interaction radius:      {:3.2f}\n'.format(CONFIG.INTERACTION_RADIUS))
        output.write('      dt:                      {:3.2f}\n'.format(CONFIG.DT))    
        output.write('      dt-beyond-i.a.-radius:   {:3.2f}\n'.format(CONFIG.DT_BEYOND_INTERA_RADIUS))     
        output.write('      mD/kT/dt !<< 1           {:3.2e}\n'.format(M_ligand / c.eta / 6e0 / math.pi /  ligand_prototype.rgyr/ CONFIG.DT *1.e10*1.e12))
        output.flush()
        
    
    
