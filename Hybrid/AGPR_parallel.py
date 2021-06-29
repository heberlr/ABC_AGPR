import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import pickle
import lhsmdu

from mpi4py import MPI
import subprocess
import os

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def CritConv(samples):
  samples=np.array(samples)
  Nsamples_a = 0.1*samples.shape[0]
  Nsamples_b = 0.5*samples.shape[0]
  Nsamples_a = int(Nsamples_a)
  Nsamples_b = int(Nsamples_b)
  Samples_a = samples[0:Nsamples_a]
  Samples_b = samples[-Nsamples_b:]
  return np.abs(np.mean(Samples_a) - np.mean(Samples_b))/np.mean(samples)

def CreateVert(Npar):
  vertices = []
  for i in range(0,2**Npar):
    temp = list(bin(i))
    temp.remove("b")
    temp.remove("0")
    while len(temp) < Npar:
        temp.insert(0,0)
    vertices.insert(0,temp)
  return np.array(vertices,dtype='f')

def CreateMesh(Gridsize, x0, xf):
  k = np.shape(x0)[0]
  M = np.zeros(shape=(Gridsize**k, k))
  d_x = np.zeros(k)
  for i in range(k):
    d = 0
    j = 0
    d_x[i] = (xf[i] - x0[i]) / (Gridsize - 1)  # increment of the cube x dimension
    x = np.arange(x0[i], xf[i]+d_x[i], d_x[i], dtype=float)
    for v in range(Gridsize ** (k - i - 1)):
      for j in range(Gridsize):
        temp = x[j]
        for z in range(Gridsize ** i):
          M[d, i] = temp
          d = d + 1
  return M

def AssPosMesh(Mesh):
  k = Mesh.shape[0]
  n = Mesh.shape[1]
  dy = Mesh[1,0] - Mesh[0,0] 
  ElemPos = []
  ElemPosX = []
  for i in range(0,k):
    for j in range(0,k):
      NB = True
      dx = Mesh[j,:] - Mesh[i,:]
      for l in range(0,n):
      	#0.1*dy fixes possible machine operation errors
        if (dx[l] > dy+0.1*dy) or (dx[l] < 0):
          NB = False
          break
      if (NB == True):
        ElemPosX.append(j)
    if (len(ElemPosX) == 2**n):
      ElemPos.append(np.array(ElemPosX))
    ElemPosX.clear()
  return np.array((ElemPos))

def ElemSampling(Ind, ParMesh, PosNode):
  n = ParMesh.shape[1] # Number of parameters
  k = PosNode.shape[1] # Number of neighbours nodes 
  M = np.zeros(n)
  d = 0
  for i in range(1,k):
    a = ParMesh[PosNode[Ind,0],:]
    b = ParMesh[PosNode[Ind,i],:]
    c = b-a
    Count = 0
    for l in range(0,n):
      if (c[l] != 0):
        Count = Count +1
    if (Count == 1):
      M[d] = ParMesh[PosNode[Ind,0],d] + np.amax(c)*np.random.uniform(0,1)
      d = d+1
  return M  

def AdapGP(Model, NsamplesInitial, LowLimit, UpperLimit, NumQOI, max_iterations=1, tol=10**-3, NSampCov=40):
  Npar = UpperLimit.shape[0]
  if rank == 0:
      samples = lhsmdu.sample(Npar,NsamplesInitial) # Latin Hypercube Sampling of Npar variables, and NsamplesInitial samples each.
      vertices = CreateVert(Npar) # Vertices of the hypercube
      samples = np.array(samples).T

      #Concatenate samples and nodes of parametric hypercube
      samples = np.concatenate((samples, vertices), axis=0)
      
      #Mesh for the parametric space
      x_0 = np.zeros(Npar)
      x_f = np.ones(Npar)

      ParMesh = CreateMesh(NsamplesInitial+1,x_0,x_f) # Parametric space unity
      PosNode = AssPosMesh(ParMesh) # Reference of the nodes associate to elemets

      #Transform unit hypercube in real values
      for j in range(0, Npar):
        samples[:,j] = LowLimit[j] + samples[:,j]*(UpperLimit[j]-LowLimit[j])
        ParMesh[:,j] = LowLimit[j] + ParMesh[:,j]*(UpperLimit[j]-LowLimit[j])
     
      #Output of the model
      OutputModel = np.zeros( (samples.shape[0],NumQOI) )
      QOI_replicas = np.zeros((size-1, NumQOI))
      for i in range(0, samples.shape[0]):     
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
      
      # Obtain prectidion for the mesh points
      GPstd = np.zeros((NsamplesInitial+1)**Npar)
      StdElem = np.zeros((NsamplesInitial)**Npar)
      for i in range(0, GPstd.shape[0]):
        OutGP = gpr.predict([np.array(ParMesh[i,:])], return_std=True)
        GPstd[i] = np.squeeze(OutGP[1])
      
      # Calculating variance in element of the mesh  
      for i in range(0, StdElem.shape[0]):
        StdElem[i] = 0
        for j in range(0, 2**Npar):
          StdElem[i] =  GPstd[PosNode[i,j]] + StdElem[i]      
        StdElem[i] = (1/(2**Npar))*StdElem[i]
      
      #Calculating difference between StdMax amd StdMean
      MeanStd = np.mean(StdElem)
      InitialSamples = samples
      MaxStdDiffMean = []
      MaxStdDiffMean.append(StdElem.max() - MeanStd)
      print("Initial step...")
      # Main loop
      for k in range(0,max_iterations):    
        StdMax = np.amax(StdElem)
        Ind = np.argmax(StdElem)
        value = ElemSampling(Ind,ParMesh,PosNode)
        for rankID in range(1,size):
            comm.Send(np.append(value,0.0), dest=rankID, tag=rankID)
        for rankID in range(1,size):
            comm.Recv(QOI_replicas[rankID-1,:], source=rankID, tag=rankID+size)
        AddSolution = np.array([np.mean(QOI_replicas, axis=0)])
        print(AddSolution.shape)
        AddSample = np.array([value])
        samples = np.concatenate((samples, AddSample), axis=0)
        OutputModel = np.concatenate((OutputModel, AddSolution), axis=0)
        
        #Generate the new fit
        gpr.fit(samples,OutputModel)
        
        # Obtaing prectidion for mesh points
        for i in range(0, GPstd.shape[0]):
          OutGP = gpr.predict([np.array(ParMesh[i,:])], return_std=True)
          GPstd[i] = np.squeeze(OutGP[1])
        
        # Calculating variance in element of the mesh  
        for i in range(0, StdElem.shape[0]):
          StdElem[i] = 0
          for j in range(0, 2**Npar):
            StdElem[i] =  GPstd[PosNode[i,j]] + StdElem[i]  
          StdElem[i] = (1/(2**Npar))*StdElem[i]

        #Calculating difference between StdMax amd StdMean
        MeanStd = np.mean(StdElem)
        MaxStdDiffMean.append(StdElem.max() - MeanStd)
        
        print("Number of new samples: "+str(k+1)+" Crit_Diff: "+str(MaxStdDiffMean[-1]))  
        #Stop criterion
        if (k >= NSampCov) and (CritConv(MaxStdDiffMean[-NSampCov:]) < tol):
          break
      
      # Finished Threads
      for rankID in range(1,size):
        comm.Send(np.append(np.zeros(Npar),1.0), dest=rankID, tag=rankID)
        
      #Finishing the metamodel
      print("==============================================================================") 
      print("Total of iterations: ", k+1) 
      print("Total of samples: ", samples.shape[0]) 
      print("Difference between Max(std) and Mean(std): ", MaxStdDiffMean[-1])
      print("==============================================================================") 
        
      # SAVE
      #GPR
      with open('model.pkl','wb') as f:
        pickle.dump(gpr,f)
      #Samples
      samples = np.matrix(samples)
      with open('samples.txt','wb') as f:
        for line in samples:
            np.savetxt(f, line, fmt='%.8e')
      #Difference of std
      np.savetxt('diffenceSTD.txt', MaxStdDiffMean, fmt='%.8e')
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
    
def Reta(Par,NumQOI):
  x = np.arange(0,1,0.01) 
  return Par[0]*x + Par[1]
  
UP = np.array([0.5,0.5])
DW = np.array([0.0,0.0])
#Training the Gaussian Process Regression (uncomment to train the GPR)
AdapGP(Calling_modelOpenMP,20, DW, UP, NumQOI=66, tol = 10e-3)

#TEST
#AdapGP(Reta,20, DW, UP, NumQOI=100, tol = 1.0)



