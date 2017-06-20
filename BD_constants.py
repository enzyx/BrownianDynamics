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
from math import pi

# Thermodynamics
k           = 1.38064850e-23       # [J/K]
eta         = 8.9e-4               # [Js/m^2]
k_per_eta   = 1.38064850/8.9 * 0.1 # [A^3/K/ps]
N_A         = 6.022e23             # [1/mol]

# Coulombic Force
e           = 1.60217662e-19 # [C]
eps0        = 8.85418781e-12 # [F/m]
ee_per_eps0 = 2.89915910e-27
epsWAT      = 82.
deg_to_rad  = pi /180.0
deg_to_rad  = pi /180.0
