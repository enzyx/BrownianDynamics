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
    
