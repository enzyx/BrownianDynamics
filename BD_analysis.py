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

def histogram(data, N_bins):
    """
    Returns the normalized histogram of data.
    The data range min(data) to max(data) is divided into N_bins. 
    """
    data_min = min(data)
    data_max = max(data)
    len_data = len(data)
    d         =  float(data_max - data_min ) / N_bins
    hist      =  numpy.zeros([N_bins,2], float)
    # sort data into histogram bins:
    for i in range(0,len(data)):
        index       = int( (data[i] - data_min) / d )
        # max(data) entry shall be included in the last bin with index N_bins-1:
        if index==N_bins:
            index = N_bins - 1
        hist[index,1] += 1
    #Assign the bin positions and normalize:
    for i in range(0, N_bins):
        hist[i,0]  = data_min + i * d
        hist[i,1] /= len_data
           
    return hist


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
