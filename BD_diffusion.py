import BD_constants as c
import math

def sigma_relative_trans(rgyr1, rgyr2, CONFIG):
    """
    Returns sigma corresponding to the diffusion coefficient
    """
    return  math.sqrt(2 * ( diffusion_trans(CONFIG.T, rgyr1) + diffusion_trans(CONFIG.T, rgyr2) ) * CONFIG.DT)

def sigma_rot(rgyr, CONFIG):
    """
    Returns sigma corresponding to the diffusion coefficient
    """
    return  math.sqrt(2 * diffusion_rot(CONFIG.T, rgyr) * CONFIG.DT)

def diffusion_trans(T, rgyr):
    """
    Returns the translational diffusion coefficient according to Stoke's law.
    """
    return c.k_per_eta * T / 6e0 / math.pi /  rgyr       # [A^2/ps]

def diffusion_rot(T, rgyr):
    """
    Returns the rotational diffusion coefficient according to Stoke's law.
    """
    return c.k_per_eta * T / 8e0 / math.pi /  rgyr**3

