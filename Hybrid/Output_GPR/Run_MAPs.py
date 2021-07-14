#! /usr/bin/env python3
#
import numpy as np
from mpi4py import MPI
import subprocess
from scipy.stats import multivariate_normal
import os

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
  
def Replicas(Model, theta, FILE='OutModel.dat'):
  Npar = theta.shape[0] # Number of parameters
  Nqoi = 66 # Number of quantity of interest
  if rank == 0:
    # Send and Receive data (MPI)
    QOI = np.zeros((size-1, Nqoi))
    for rankID in range(1,size):
        comm.Send(np.append(theta,0.0), dest=rankID, tag=rankID)
    for rankID in range(1,size):
        comm.Recv(QOI[rankID-1,:], source=rankID, tag=rankID+size)
    
    # Write in file
    np.savetxt(FILE, QOI,delimiter='\t', fmt='%e')
   
    # Finished Threads
    for rankID in range(1,size):
        comm.Send(np.append(theta,1.0), dest=rankID, tag=rankID)
  else:
      Par = np.zeros(Npar+1, dtype='d')   
      while (Par[-1]==0.0):
        comm.Recv(Par, source=0,tag=rank)
        if (Par[-1] == 1.0): break
        OUT = Model(Par[:-1],Nqoi)
        comm.Send(OUT, dest=0,tag=rank+size)
  
def Calling_modelOpenMP(parameter,NumQOI):
    function_call = ['./build/main.exe', '{}'.format(rank)] #change of .exe to .out
    for ind in range(0,parameter.shape[0]):
        function_call.append('{}'.format(parameter[ind]))
    cache = subprocess.run(function_call,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(cache.stdout)
    # print(cache.stderr)
    # print(cache.returncode)
    #print(OUT)
    OUT = np.fromstring(cache.stdout, dtype='d', sep=' ')
    time = OUT[::3]
    Live = OUT[1::3]
    Dead = OUT[2::3]
    QOI = np.concatenate((Live, Dead), axis=None)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)
    if (QOI.shape[0] != NumQOI):
        print("Model output error! incompatibility of QoIs!")
        os._exit(0)
    return QOI  
    
theta_star1 = np.array([0.2098,0.1632], dtype='d') # 100_MCMC_GPR
theta_star2 = np.array([0.2432,0.1481], dtype='d') # 100_SMC_GPR
theta_star3 = np.array([0.2614,0.1554], dtype='d') # 1000_MCMC_GPR
theta_star4 = np.array([0.2127,0.1686], dtype='d') # 1000_SMC_GPR

Replicas(Calling_modelOpenMP,theta_star1,"MAPs/OutModel_100_MCMC_GPR.dat")
Replicas(Calling_modelOpenMP,theta_star2,"MAPs/OutModel_100_SMC_GPR.dat")
Replicas(Calling_modelOpenMP,theta_star3,"MAPs/OutModel_1000_MCMC_GPR.dat")
Replicas(Calling_modelOpenMP,theta_star4,"MAPs/OutModel_1000_SMC_GPR.dat")
