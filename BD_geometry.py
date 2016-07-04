import math
import numpy
from math import cos, sin

def translate(molecule, dR):
    """
    translate all atom coordinates in R 
    """
    for i in range(molecule.N_atoms):
        molecule.R[i] += dR
    molecule.center += dR


def rotate(molecule, dO, rotate_around_center = False):
    """
    rotate all atom coordinates in R around the origin
    """
    sinalpha       = sin(dO[0])
    cosalpha       = cos(dO[0])
    rot_alpha      = numpy.matrix([[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]])
    
    rot_alpha[0,0] = cosalpha
    rot_alpha[0,1] = sinalpha
    rot_alpha[1,0] = -sinalpha
    rot_alpha[1,1] = cosalpha    

    sinbeta        = sin(dO[1])
    cosbeta        = cos(dO[1])
    rot_beta   = numpy.matrix([[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]])
    
    rot_beta[1,1]  = cosbeta
    rot_beta[1,2]  = sinbeta
    rot_beta[2,1]  = -sinbeta
    rot_beta[2,2]  = cosbeta

    singamma       = sin(dO[2])
    cosgamma       = cos(dO[2])
    rot_gamma  = numpy.matrix([[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]])   
    
    rot_gamma[0,0] = cosgamma
    rot_gamma[0,1] = singamma
    rot_gamma[1,0] = -singamma
    rot_gamma[1,1] = cosgamma

    for i in range(molecule.N_atoms):
        if rotate_around_center == True:
            R_matrix = numpy.matrix([[molecule.R[i,0]] - molecule.center[0],[molecule.R[i,1]] - molecule.center[1], [molecule.R[i,2] - molecule.center[2]]]) 
            R_matrix = rot_gamma * rot_beta * rot_alpha * R_matrix
            for d in range(3):
                molecule.R[i,d] = R_matrix[d,0] + molecule.center[d]
        else: 
            R_matrix = numpy.matrix([[molecule.R[i,0]], [molecule.R[i,1]], [molecule.R[i,2]]]) 
            R_matrix = rot_gamma * rot_beta * rot_alpha * R_matrix
            for d in range(3):
                molecule.R[i,d] = R_matrix[d,0]
    
    if rotate_around_center == False:
        R_matrix = numpy.matrix([[molecule.center[0]],[molecule.center[1]],[molecule.center[2]]]) 
        R_matrix = rot_gamma * rot_beta * rot_alpha * R_matrix   
        for d in range(3):
            molecule.center[d] = R_matrix[d,0]   

def center(R, indices = 0):
    """
    returns the geometric center of all atom positions in R
    if a indices list is given, returns the geometric center of 
    the atoms corresponding to the indices
    """
    center = numpy.zeros([3])
    if indices == 0 or indices == [0]:
        for index in range(len(R)):
            center += R[index]
        return center / len(R)
    else:
        for index in indices:
            center += R[index]
        return center / len(indices)        


def rgyr(molecule):
    """
    calculates the radius of gyration of a molecule
    """
    return   numpy.sqrt( numpy.sum( (molecule.R - molecule.center)**2 ) / molecule.N_atoms )

def rmax(molecule):
    """
    returns the largest atom distance from the center of the molecule 
    """
    V = numpy.zeros([molecule.N_atoms])
    for i in range(molecule.N_atoms):
        V[i] = (molecule.R[i,0] - molecule.center[0])**2 + (molecule.R[i,1] - molecule.center[1])**2 + (molecule.R[i,2] - molecule.center[2])**2
    return numpy.sqrt( numpy.max(V) )
          
 
# Cartesian -> Spherical
def r(R):
    return math.sqrt(R[0]**2 + R[1]**2 +R[2]**2) 

def theta(R):
    return math.acos( R[2] / math.sqrt(R[0]**2 + R[1]**2 +R[2]**2) )

def phi(R):
    if   R[0] > 0.0:                      
        return math.atan( R[1] / R[0] )
    elif R[0] == 0.0:
        return numpy.sign(R[1]) * 0.5 * math.pi 
    elif R[0] < 0.0 and R[1] >= 0.0:
        return math.atan( R[1] / R[0] ) + math.pi
    elif R[0] < 0.0 and R[1] < 0.0:
        return math.atan( R[1] / R[0] ) - math.pi

# def S_spherical(R):
#     S    = numpy.zeros([3], dtype = 'float')
#     S[0] = r(R)
#     S[1] = theta(R)
#     S[2] = phi(R)
#     return S
# 
# # Spherical -> Cartesian    
# def x(S):
#     return S[0] * math.sin(S[1]) * math.cos(S[2])
# 
# def y(S):
#     return S[0] * math.sin(S[1]) * math.sin(S[2])
# 
# def z(S): 
#     return S[0] * math.cos(S[1])
# 
# def R_cartesian(S):
#     R    = numpy.zeros([3], dtype = 'float')
#     R[0] = x(S)
#     R[1] = y(S)
#     R[2] = z(S)
#     return R   
#===============================================================================
    
