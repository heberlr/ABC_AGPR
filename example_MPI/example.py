import numpy as np
import pickle #library to store the model trained

# Changing repository
import sys
sys.path.append('../src')
# Loading function to ABC calibration
from ABC_parallel import ABC_SMC,ABC_MCMC
from AGPR_parallel import AdapGP


# MPI initializing
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def GenerateData():
  x = np.linspace(0,1,100)
  Par = np.array([2.0,0.5])
  return Par[0]*x + Par[1]
  
def Reta(Par,Nqoi):
  x = np.linspace(0,1,Nqoi)
  return Par[0]*x + Par[1] + np.random.normal( 0.0,0.01) # Noise from normal dist.

if __name__ == '__main__':
  # Observational data
  Data = GenerateData()
  # Boundary of parameter space
  UpperLimit = np.array([3.0,1.0])
  LowLimit = np.array([1.0,0.0])
  # Tolerance vector
  epsilon = np.array([1.0,0.4,0.2])

  #ABC_MCMC(Reta,Data,LowLimit,UpperLimit,Nrep=10,tol=epsilon[-1],NumAccept=100)
  #ABC_SMC(Reta,Data,LowLimit,UpperLimit,Nrep=10,tol=epsilon,NumAccept=100)
  
  #Training the Gaussian Process Regression (uncomment to train the GPR)
  NsamplesInitial=20
  AdapGP(Reta,NsamplesInitial, LowLimit, UpperLimit, Data.size, NumRep =10, tol = 1e-3)