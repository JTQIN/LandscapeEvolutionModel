import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import warnings
import os
import FlowRouting
import shutil
import timeit
warnings.filterwarnings('error')

output_folder = 'Z:\\ACTIVE\\MiladHooshyar\\LEM\\Test\\'

D = 0.031461158
K = 0.000555381
U = 0.0001
m = 0.324
dX = 5
n = 1
dt_min = 116

T_save = 5000
min_Err = 0.01
max_Err = 20.
dZ1 = 0.05
dZ2 = 0.1
min_T_conv = 10000

Nx = 301
Ny = 301
H = 1.
    
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)

os.makedirs(output_folder)
X = np.arange(0 , Nx * dX , dX)
Y = np.arange(0 , Ny * dX, dX)
Con_Data = []
## Initial condition
H_min = 1000.
dH = 2. * H / Ny
Z_0 = np.zeros((Ny , Nx))
Z_1 = np.zeros((Ny , Nx))
for i in range(0 , Ny):
    if i <= Ny/2:
        Z_0[i , :] = H_min + i * dH
    else:
        Z_0[i , :] = H_min + 2 * H - i * dH
np.random.seed(seed=700)
Z_0 = Z_0 + np.random.normal(0, 0.01, [Ny,Nx])
Z_0 = np.where(Z_0 < H_min , H_min , Z_0)

## Simulation
j = 0
t = 0.
Err = 10.
ct_save = 0
dt = 0
T_conv = 0
while  T_conv < min_T_conv or t <= 5 * 10**6:
    t = t +  dt
    start = timeit.default_timer()
    Z_0 , A, fdir , pit = FlowRouting.flow_routing(Z_0 , dX , 'NoFill')
    stop = timeit.default_timer()
    
    A = A * dX**2
    lap , slp = FlowRouting.slope_lap(Z_0 , dX)

    dZ = (D * lap - K * A ** m * slp * (1 - pit) + U)
    max_dZ = np.max(np.abs(dZ))
    if t <= 10000:
        dt = max(min(int(dZ1 / max_dZ) , 1000) , int(0.5 * dt_min))
    else:
        dt = max(min(int(dZ2 / max_dZ) , 1000) , int(dt_min))
    
    Z_1 = dZ * dt + Z_0
    Z_1 = np.where(Z_1 < H_min + 0.01 , H_min + 0.01 , Z_1)
    Z_1[0, :] = H_min
    Z_1[Ny - 1, :] = H_min
    Con_Data.append([0 , t , np.sum(np.abs(Z_1 - Z_0))/dt , np.sum((Z_1 - Z_0))/dt , round(np.max(np.abs(Z_1 - Z_0))/dt * 1000 , 2) , round(stop - start , 2), int(np.sum(fdir==-1))])
    Err = abs(round(np.sum((Z_1 - Z_0))/dt , 3))
    if Err <= min_Err:
        T_conv = T_conv + dt
    else:
        T_conv = 0
    if Err >= max_Err:
        listNoConv.append([tr , 0])
        break
    if t >= 3 * 10**6:
        listNoConv.append([tr , 1])
        break

    print ' year= ', int(t) , ', dt= ', int(dt) , ', total ele. change= ', Err , ', run time= ' , round(stop - start , 2),'s ' , ', num of sinks= ' , int(np.sum(fdir==-1))
    if t <= (ct_save * T_save) and (t + dt) > (ct_save * T_save):
        out_array = output_folder + '\\DEM_' + str(int(ct_save)) +'.npy'
        np.save(out_array, Z_1)
        
        out_DEM = output_folder + '\\DEM_' + str(int(ct_save)) +'.tif'
        arcpy.NumPyArrayToRaster(Z_1).save(out_DEM)
        
        Data = np.asarray(Con_Data)
        ct_save = ct_save + 1
    Z_0 = np.copy(Z_1)
    j = j + 1

