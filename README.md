# extragas_grid

Building extragas models in parallel. 


# Outline: 

Model bulding:  C++ codes (extragas)       

We wrap the C++ codes into a shared library. 

Use Python codes to call the shared library and do the parallel computing. 


# Steps:

(1) Go to directory "extragas_grid/c++".   

(2) Compile the C++ codes 

$g++ -mtune=native  -g -Wall -o /data/users/li/extragas_grid/library/n2403_model.so -shared -fPIC fitsUtils.cpp fitsutil.cpp myutilities.cpp statist.cpp numint.cpp interpol.cpp  findValue.cpp  galFunctions.cpp dynamics.cpp extragasFits.cpp extragasPlots.cpp projections.cpp iterations.cpp inout.cpp accretion.cpp hothalo.cpp extragas.cpp extragas_fit_n2403.cpp main_wrap.cpp  -w -L/data/users/li/extragas_grid/library -lcfitsio -O3 

"/data/users/li/extragas_grid/library/n2403_model.so": name and directory of the output library name 
 
"/data/users/li/extragas_grid/library": PATH for the compiler to search lcfitsio library  


the output is a shared library n2403_model.so for python codes. 


(3)Go to directory "extragas_grid/python".

(4)Edit grid_run.py 

result_path: path for storing output models fits files. 

libname:   the name and directory of the shared library generated in the previous step. 
 
(5) start the parallel computing: 

$mpirun -hostfile hostfile --map-by ppr:30:node python3 grid_run.py

hostfile is a file listing all the available cluster nodes 

ppr:N_process:node: N_process defines how many processes to run on each node. usually it should be smaller than the N_cpu on each node. 





