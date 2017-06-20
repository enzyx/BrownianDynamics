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
import parmed
import numpy
import BD_geometry as g
import sys

def load_molecule_from_pqr(filename):
    """
    -loads a molecule from a .pdb file using parmed.
    -converts the atom positions into a list of one-column matrices
    -converts the atom charges into a list of floats
    """
    
    mol = parmed.load_file(filename)
    
    R = numpy.zeros([len(mol.atoms),3], 'float')
    q = numpy.zeros([len(mol.atoms)], 'float')
    
    l = len(R)
    
    for i in range(len(R)):
        sys.stdout.write('\r   loading {} atom {:6d}/{:6d}'.format(filename, i+1, l))
        sys.stdout.flush()
        for d in range(3):
            R[i,d] = mol.coordinates[i,d]
        q[i] = mol.atoms[i].charge
    
    sys.stdout.write('\n')
    sys.stdout.flush()
    return R, q
    
def write_molecules_to_pqr(molecule, filename, original_pqr_filename, CONFIG):
    mol = parmed.load_file(original_pqr_filename)
    for i in range(molecule.N_atoms):
        mol.atoms[i].xx = molecule.R[i,0]
        mol.atoms[i].xy = molecule.R[i,1]
        mol.atoms[i].xz = molecule.R[i,2]
    
    parmed.write_PDB(mol, CONFIG.OUTPUT_DIRECTORY+'/'+filename+'.pqr')
    
