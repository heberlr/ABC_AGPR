# mpiexec -n 4 py -m mpi4py RunModel.py

import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import pickle
import lhsmdu
import sys

from mpi4py import MPI
import subprocess
import os

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def GPR(Model, Nsamples, LowLimit, UpperLimit, NumQOI):
  Npar = UpperLimit.shape[0]
  if rank == 0:
      samples = np.random.uniform(0,1,Nsamples*Npar)
      samples = np.reshape(samples, (Nsamples,Npar))

      #Transform unit hypercube in real values
      for j in range(0, Npar):
        samples[:,j] = LowLimit[j] + samples[:,j]*(UpperLimit[j]-LowLimit[j])
     
      #Output of the model
      OutputModel = np.zeros( (samples.shape[0],NumQOI) )
      QOI_replicas = np.zeros((size-1, NumQOI))
      for i in range(0, samples.shape[0]): 
        print("Sample: "+str(i)+" Par: "+str(samples[i,:]))
        for rankID in range(1,size):
            comm.Send(np.append(samples[i,:],0.0), dest=rankID, tag=rankID)
        for rankID in range(1,size):
            comm.Recv(QOI_replicas[rankID-1,:], source=rankID, tag=rankID+size)
        OutputModel[i,:] = np.mean(QOI_replicas, axis=0)   

      #Set initial hyperparemeter of the GPR
      hyp = np.zeros(Npar)
      for j in range(0, Npar):
        hyp[j] = 0.1*(UpperLimit[j]-LowLimit[j])
      
      #Set kernel
      kernel = RBF(hyp)
      
      #Set GPR options
      gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=50, alpha=10**-5)
      
      #Generate the initial fit
      gpr.fit(samples,OutputModel)
      
      # Finished Threads
      for rankID in range(1,size):
        comm.Send(np.append(np.zeros(Npar),1.0), dest=rankID, tag=rankID)
        
      #Finishing the metamodel
      print("==============================================================================") 
      print("Total of samples: ", samples.shape[0])
      print("==============================================================================") 
        
      # SAVE
      #GPR
      fileModel = 'GPR_'+"%03i"%Nsamples+'/modelGPR.pkl'
      with open(fileModel,'wb') as f:
        pickle.dump(gpr,f)
      #Samples
      samples = np.matrix(samples)
      fileSamples = 'GPR_'+"%03i"%Nsamples+'/samplesGPR.txt'
      with open(fileSamples,'wb') as f:
        for line in samples:
            np.savetxt(f, line, fmt='%.8e')
  else:
      Par = np.zeros(Npar+1, dtype='d')   
      while (Par[-1]==0.0):
        comm.Recv(Par, source=0,tag=rank)
        if (Par[-1] == 1.0): break
        OUT = Model(Par[:-1],NumQOI)
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
  
UP = np.array([0.5,0.5])
DW = np.array([0.0,0.0])
#Training the Gaussian Process Regression (uncomment to train the GPR)
GPR(Calling_modelOpenMP,int(sys.argv[1]), DW, UP, NumQOI=66)

#TEST
#AdapGP(Reta,20, DW, UP, NumQOI=100, tol = 1.0)



