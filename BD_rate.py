#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
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
