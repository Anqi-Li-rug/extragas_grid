import ctypes
import mpi4py
import sys
import time
import os
import numpy as np
from schwimmbad import MPIPool
from astropy.io import fits
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import emcee
import corner

#path to store the output model cubes. Note it has to be the absolute path (relative to the whole clusters) since we will be running the code on several machines.
result_path='/net/cuby3.apertif/data/users/li/model_cubes/'
#path and name of the c++ library
libname="/net/cuby3.apertif/data/users/li/extragas_grid/library/n2403_model.so"

#setting the steps and ranges of the free-parameter-grid;  

# ejection velocity, [km/s]
min_v_k=40.
max_v_k=110.
step_v_k=10.

#ionisation fraction
min_f_ion=0.0
max_f_ion=1.0
step_f_ion=0.2

# condensataion rate alpha,   [1/Gyr]
min_alpha=0.0
max_alpha=6.0
step_alpha=0.6


#Call c++ code from python;
#library loading
c_lib = ctypes.cdll.LoadLibrary(libname)
#setting types for the arguments
c_lib.grid.restype = ctypes.c_double
c_lib.argtypes = [ctypes.c_float,ctypes.c_float,ctypes.c_float, ctypes.c_char_p]



def grid(task):
    dem=True
    while dem==True:
        try:
            v_k,f,accre,fitsname=task
            if os.path.exists(fitsname):
                print("this model has already been generated")
                return 0
            print(v_k,f,accre,fitsname)
            z=c_lib.grid(ctypes.c_float(v_k), ctypes.c_float(f),ctypes.c_float(accre),fitsname) 
            if z<0.0000001:
                print(v_k,f,accre)
                print("something went wrong with this param!")
                dem=True
            else:
                dem=False
        except:
            dem=True
    return z
def main(pool):
    #here we calculate the parameter values for each grid point
    v_k=np.arange(min_v_k,max_v_k+step_v_k,step_v_k)
    f_ion=np.arange(min_f_ion,max_f_ion+step_f_ion,step_f_ion)
    alpha=np.arange(min_alpha,max_alpha+step_alpha,step_alpha)

    # wrap the free parameters and the output model file names into a format readable by C++ function
    name_list=[]
    f_ion_list=[]
    v_k_list=[]
    alpha_list=[]
    for i in np.arange(len(v_k)):
        for j in np.arange(len(f_ion)):
            for k in np.arange(len(alpha)):
                v_k_list.append(v_k[i])
                f_ion_list.append(f_ion[j])
                alpha_list.append(alpha[k])
                name=result_path+'v'+str('{:.0f}').format(v_k[i])+'_fion'+str('{:.1f}').format(f_ion[j])+'_alpha'+str('{:.2f}').format(alpha[k])+'.fits'
                name_encode=name.encode('utf-8')
                name_list.append(name_encode)

    #arguments passing to the c++ function
    tasks = list(zip(v_k_list,f_ion_list,alpha_list,name_list))
    #parallel modelling
    results = pool.map(grid, tasks)
    pool.close()

#main   pool
if __name__ == "__main__":
    import sys
    from schwimmbad import MPIPool
    from mpi4py import MPI
    pool=MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    main(pool)
