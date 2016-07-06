#!/usr/bin/python2
# -*- coding: utf-8 -*-
import numpy
import sys
import argparse

import BD_molecules as m
import BD_configparser as config
import BD_diffusion as diffusion
import BD_constants as c

# 1. INPUT 
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Calculate the k_on rate into the encounter complex.')
parser.add_argument('-c','--config-file', type=str, dest="CONFIGFILE",
                    required=False, default = 'BD.conf',
                    help="Configfile of Brownian Dynamics Simulation")

args    = parser.parse_args()
CONFIG  = config.Configuration(args.CONFIGFILE)

# 2. Load data 
# -----------------------------------------------------------------------------
# beta_BD from simulations folder
try:
    beta_BD = numpy.loadtxt(CONFIG.OUTPUT_DIRECTORY + "/end_states.cum")[-1]
except:
    print "Error: Could not find file " + CONFIG.OUTPUT_DIRECTORY + "/end_states.cum"
    print "       Maybe the simulation has not yet been run?"
    sys.exit(-1)

receptor = m.Receptor(CONFIG)
ligand_prototypes = []
for ligand_pqr in CONFIG.LIGAND_PQRS:
    ligand_prototypes.append(m.Ligand(ligand_pqr, None))

ligand_rgyr_mean = numpy.array([x.rgyr for x in ligand_prototypes]).mean()

D1 = diffusion.diffusion_trans(CONFIG.T, receptor.rgyr)
D2 = diffusion.diffusion_trans(CONFIG.T, ligand_rgyr_mean)
# Relative Diffusion constant for two diffusing spheres (see Northrup)
D = D1 + D2                       # [A^2/ps]
b = CONFIG.STARTING_RADIUS        # [A]
q = CONFIG.MAXIMUM_RADIUS         # [A]

# 3. Calculate Rate
# -----------------------------------------------------------------------------
def calculateAssociationRate():
    """
    Formula taken from
    S. H. Northrup, J. Chem. Phys. 80 (4), 15 February 1984
    """
    beta_inf = beta_BD / (1. - (1.-beta_BD) * b/q)
    return 4. * numpy.pi * D * b * beta_inf * c.N_A * 1e-15  # 1/(Ms)

k_on = calculateAssociationRate()
    
print "k_on(inf->ES) = {:e} [s⁻¹M⁻¹]".format(k_on)