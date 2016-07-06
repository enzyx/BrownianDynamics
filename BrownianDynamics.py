#!/usr/bin/python2
import numpy
import sys
import argparse
from multiprocessing import Pool, Manager
import time
import os
import cPickle as pickle
import random
import copy

import BD_geometry as g
import BD_forces as f
import BD_propagate as p
import BD_molecules as m
import BD_pqr as pqr
import BD_configparser as config

# 1. INPUT 
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Brownian Dynamics Simulation.')
parser.add_argument('-c','--config-file', type=str, dest="CONFIGFILE",
                    required=False, default = 'BD.conf',
                    help="Configfile for Brownian Dynamics Simulation")

args    = parser.parse_args()
CONFIG  = config.Configuration(args.CONFIGFILE)


# 2. PREPARATION
# -----------------------------------------------------------------------------
sys.stdout.write('\033[1mBrownian Dynamics\033[0m\n')

# Directory
if os.path.exists(CONFIG.OUTPUT_DIRECTORY):
    sys.stdout.write(' Error: output directory exists\n')
    sys.stdout.flush()
    sys.exit()
else:
    os.mkdir(CONFIG.OUTPUT_DIRECTORY) 


# load from pickle
if (not CONFIG.LOAD_MOLS_FROM_PICKLE.lower() in ['false', 'no']) and \
   (os.path.isfile(CONFIG.LOAD_MOLS_FROM_PICKLE)):
    sys.stdout.write(' Loading ligand, receptor and grid from pickle...\n')
    sys.stdout.flush()
    pickle_file = open(CONFIG.LOAD_MOLS_FROM_PICKLE, 'r')
    receptor, grid, ligand_prototypes = pickle.load(pickle_file)
    pickle_file.close()
# set up
else:
    if not os.path.isfile(CONFIG.LOAD_MOLS_FROM_PICKLE):
        sys.stdout.write(' No pickle file found...\n')
    sys.stdout.write(' Setting up receptor...\n')
    sys.stdout.flush()
    receptor            = m.Receptor(CONFIG)
    
    sys.stdout.write(' Setting up grid...\n')
    sys.stdout.flush()
    grid                = m.Grid(receptor, CONFIG)
    
    sys.stdout.write(' Setting up ligand...\n')
    sys.stdout.flush()
    
    ligand_prototypes = []
    for ligand_pqr in CONFIG.LIGAND_PQRS:
        ligand_prototypes.append(m.Ligand(ligand_pqr, grid))
    
    if not (CONFIG.SAVE_MOLS_TO_PICKLE.lower() in ['false', 'no']):  
        pickle_file = open(CONFIG.SAVE_MOLS_TO_PICKLE, 'w')
        pickle.dump([receptor, grid, ligand_prototypes], pickle_file)
        pickle_file.close()


pqr.write_molecules_to_pqr(receptor, 'receptor_centered', CONFIG.RECEPTOR_PQR, CONFIG)    
propagator = p.Propagator(ligand_prototypes[0], receptor, CONFIG)
config.checkInputConsistency(ligand_prototypes[0], receptor, grid, propagator, CONFIG)
config.printSimulationInfo(ligand_prototypes[0], receptor, grid, propagator, CONFIG)


# 3. PROPAGATE TRAJECTORIES
# -----------------------------------------------------------------------------
sys.stdout.write(' Propagation...\n')
sys.stdout.flush()

def trajectory(args):
    i, queue = args
    # Preparation -------------------------------------------------------------
    # initiate a separate instance of Random for every trajectory,
    # so that in the case of threading every thread gets different random numbers
    ran = random.Random()
    ran.jumpahead(CONFIG.NUMBER_TRAJECTORIES)
    
    if CONFIG.RANDOM_START_POSITIONS == True:
        # create an own copy of the propagated object for every trajectory run
        ligand = copy.deepcopy(ligand_prototypes[0])
        # move it to a random position on a sphere of radius STARTING_RADIUS
        # with random orientation and update grid indices
        propagator.starting_position(ligand, CONFIG, ran)
    else:
        # create an own copy of the propagated object for every trajectory run
        ligand = copy.deepcopy(ligand_prototypes[i%len(CONFIG.LIGAND_PQRS)])
        #print i%len(CONFIG.LIGAND_PQRS), (CONFIG.LIGAND_PQRS), '\n'
        # translate it with the receptor to the system origin
        g.translate(ligand, - receptor.center_original)
    PHI_start   = g.phi(ligand.center)
    THETA_start = g.theta(ligand.center)         

        

    
    # check initial state of the system
    state_tmp = propagator.state(ligand, receptor, CONFIG)
    if not state_tmp ==2:
        print '\nError: Ligand starting position not between rmax and receptor'

    # Check initial ligand position for collision
    # calculate initial forces
    collision       = f.collision(ligand, receptor, grid)
    if collision:
	print '\nWARNING: starting position collides... '
    while collision:
        print 'adjusting\n'
        g.translate(ligand, 0.1*ligand.center)
    	collision       = f.collision(ligand, receptor, grid)
    
    #pqr.write_molecules_to_pqr(ligand, 'starting'+str(i), CONFIG.LIGAND_PQR, CONFIG)
    ligand.updateGridIndices(grid)
    interacting     = f.interacting(ligand, CONFIG)
    F, T_lig, T_rec = f.calcForceAndTorque(ligand, receptor, grid, CONFIG)
    
    if CONFIG.SAVE_START_COORDS:
        pqr.write_molecules_to_pqr(ligand, 'start_{:06d}'.format(i), CONFIG.LIGAND_PQRS[0], CONFIG)
    
    # Propagation Loop --------------------------------------------------------
    collision_counter = 0

    while state_tmp == 2:
    
        ligand.savePreviousPosition()
        propagator.propagate_trans(ligand, ran, interacting, F)
        propagator.propagate_rot_ligand(ligand, ran, interacting, T_lig)
        propagator.propagate_rot_receptor(ligand, ran, interacting, T_rec)
        ligand.updateGridIndices(grid)

        collision = f.collision(ligand, receptor, grid)
        if collision:
            ligand.resetToPreviousPosition()
            collision_counter += 1
            if collision_counter > 10000:
                pqr.write_molecules_to_pqr(ligand, 'collision_{:06d}'.format(i), CONFIG.LIGAND_PQRS[0], CONFIG)
                return 3, 0, 0, 0, 0 
            continue
        
        collision_counter = 0
        
        interacting = f.interacting(ligand, CONFIG)
        if interacting: 
            F, T_lig, T_rec = f.calcForceAndTorque(ligand, receptor, grid, CONFIG)
            

        state_tmp = propagator.state(ligand, receptor, CONFIG)
        if state_tmp == 1:
            pqr.write_molecules_to_pqr(ligand, 'target_{:06d}'.format(i), CONFIG.LIGAND_PQRS[0], CONFIG)

    # -------------------------------------------------------------------------
    PHI_end   = g.phi(ligand.center)
    THETA_end = g.theta(ligand.center) 
    
    if not queue == None:
        queue.put(state_tmp)
    return state_tmp, PHI_start, THETA_start, PHI_end, THETA_end


# SERIAL
if CONFIG.NUMBER_THREADS == 1:
    end_states = []
    for i in range(CONFIG.NUMBER_TRAJECTORIES):
        args = (i, None)
        end_states.append(trajectory(args))
        
        # Monitor progress
        if i%10 == 0:
            progress_file = open(CONFIG.OUTPUT_DIRECTORY+'/progress', 'w')
            progress_file.write('Trajectory:   {:7d}/{:7d}\n'.format(queue.qsize(), CONFIG.NUMBER_TRAJECTORIES))
            progress_file.close() 
                
# PARALLEL (THREADS)
elif CONFIG.NUMBER_THREADS > 1:
    pool    = Pool(processes = CONFIG.NUMBER_THREADS)
    manager = Manager()
    queue   = manager.Queue()
    args = [(i, queue) for i in range(CONFIG.NUMBER_TRAJECTORIES ) ]
    result = pool.map_async(trajectory, args)
    
    # Monitor progress
    while True:
        if result.ready():
            break
        else:
            progress_file = open(CONFIG.OUTPUT_DIRECTORY+'/progress', 'w')
            progress_file.write('Trajectory:   {:7d}/{:7d}\n'.format(queue.qsize(), CONFIG.NUMBER_TRAJECTORIES))
            progress_file.close() 
            time.sleep(10)

    end_states = result.get()

# 4. ANALYSIS AND OUTPUT
# -----------------------------------------------------------------------------
sys.stdout.write('\r  Trajectory:   {:6d}/{:6d}\n'.format(CONFIG.NUMBER_TRAJECTORIES, CONFIG.NUMBER_TRAJECTORIES))
sys.stdout.flush() 

def cumulative_mean(data):
    """
    Cumulative mean of data
    """
    cumulative_mean = []
    cumulative_sum = 0.0        
    for i, value in enumerate(data):
        cumulative_sum += value
        cumulative_mean.append(cumulative_sum/float(i+1))
    return cumulative_mean

outfile = open(CONFIG.OUTPUT_DIRECTORY+'/end_states.dat', 'w')
for i in range(len(end_states)):
    outfile.write('{:d} {:+5.2e} {:+5.2e} {:+5.2e} {:+5.2e}\n'.format(end_states[i][0], 
                                                                end_states[i][1], 
                                                                end_states[i][2], 
                                                                end_states[i][3], 
                                                                end_states[i][4]))
outfile.close()

end_states_array = []
for i in range(len(end_states)):
    if not end_states[i][0] == 3:
        end_states_array.append(end_states[i][0])
cumulative = cumulative_mean(end_states_array)
numpy.savetxt(CONFIG.OUTPUT_DIRECTORY+'/end_states.cum', cumulative)

sys.stdout.write('Fraction of trajectories arrived at target: {}\n'.format(numpy.mean(end_states_array)))
sys.stdout.flush() 
