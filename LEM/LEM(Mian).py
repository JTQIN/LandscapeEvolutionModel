import numpy as np
import FlowRouting
import timeit
import arcpy
import os
import shutil

output_folder = 'Z:\\ACTIVE\\MiladHooshyar\\LEM\\Test\\'
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.makedirs(output_folder)
    

D = 0.031461158
K = 0.000555381
U = 0.0001
m = 0.324
n = 1

dX = 5
dt = 120

Nx = 301
Ny = 301

X = np.arange(0 , Nx * dX , dX)
Y = np.arange(0 , Ny * dX, dX)

H = 1.
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
t = 0.
j = 0
while  t <= 3 * 10 ** 6:
    t = t +  dt
    start = timeit.default_timer()
    Z_0 , A, fdir , pit = FlowRouting.flow_routing(Z_0 , dX , 'NoFill')
    stop = timeit.default_timer()
    
    A = A * dX**2
    lap , slp = FlowRouting.slope_lap(Z_0 , dX)

    dZ = (D * lap - K * A ** m * slp * (1 - pit) + U)

    Z_1 = dZ * dt + Z_0
    Z_1 = np.where(Z_1 < H_min + 0.01 , H_min + 0.01 , Z_1)
    Z_1[0, :] = H_min
    Z_1[Ny - 1, :] = H_min

    print ' year= ', int(t) , ', Mean ele= ',  round(np.mean(Z_1) , 2) , ', run time= ' , round(stop - start , 2),'s ' , ', num of sinks= ' , int(np.sum(fdir==-1))

    
    if j % 5 == 0:
        out_array = output_folder + 'DEM_' + str(int(j)) +'.npy'
        np.save(out_array, Z_1)
        
        out_DEM = output_folder + 'DEM_' + str(int(j)) +'.tif'
        arcpy.NumPyArrayToRaster(Z_1).save(out_DEM)
    Z_0 = np.copy(Z_1)
    j += 1
        

