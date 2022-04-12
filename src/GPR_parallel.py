import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import pickle

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

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

def ElemSampling(Npar, LowLimit, UpperLimit):
  sample = np.zeros(Npar, dtype='d')
  for j in range(0, Npar):
    sample[j] = LowLimit[j] + np.random.uniform(0,1)*(UpperLimit[j]-LowLimit[j])
  return sample

def GPR(Model, NpartionsLHD, LowLimit, UpperLimit, NumQOI, folder, NumRep = 10, max_samples=1000):
    Npar = UpperLimit.shape[0]
    if rank == 0:
        samples = np.random.uniform(0,1,NpartionsLHD*Npar)
        samples = np.reshape(samples, (NpartionsLHD,Npar))
        #Mesh for the parametric space
        x_0 = np.zeros(Npar)
        x_f = np.ones(Npar)
        ParMesh = CreateMesh(NpartionsLHD+1,x_0,x_f) # Parametric space unity
        PosNode = AssPosMesh(ParMesh) # Reference of the nodes associate to elemets
        #Transform unit hypercube in real values
        for j in range(0, Npar):
            samples[:,j] = LowLimit[j] + samples[:,j]*(UpperLimit[j]-LowLimit[j])
            ParMesh[:,j] = LowLimit[j] + ParMesh[:,j]*(UpperLimit[j]-LowLimit[j])
        #Output of the model
        OutputModel = np.zeros( (samples.shape[0],NumQOI) )
        QOI_replicas = np.zeros((NumRep, NumQOI))
        for i in range(0, samples.shape[0]):
            for indx in range(0,NumRep):
                rankID = indx%(size-1) + 1
                comm.Send(np.append(samples[i,:],0.0), dest=rankID, tag=rankID)
            for indx in range(0,NumRep):
                rankID = indx%(size-1) + 1
                comm.Recv(QOI_replicas[indx,:], source=rankID, tag=rankID+size)
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
        GPstd = np.zeros((NpartionsLHD+1)**Npar)
        StdElem = np.zeros((NpartionsLHD)**Npar)
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
        for k in range(NpartionsLHD,max_samples):
            StdMax = np.amax(StdElem)
            Ind = np.argmax(StdElem)
            value = ElemSampling(Npar, LowLimit, UpperLimit)
            # Send and Receive data (MPI)
            for rankID in range(0,NumRep):
                rankID = indx%(size-1) + 1
                comm.Send(np.append(value,0.0), dest=rankID, tag=rankID)
            for rankID in range(0,NumRep):
                rankID = indx%(size-1) + 1
                comm.Recv(QOI_replicas[rankID-1,:], source=rankID, tag=rankID+size)
            AddSolution = np.array([np.mean(QOI_replicas, axis=0)])
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

        # Finished Threads
        for rankID in range(1,size):
            comm.Send(np.append(np.zeros(Npar),1.0), dest=rankID, tag=rankID)

        #Finishing the metamodel
        print("==============================================================================")
        print("Total of samples: ", samples.shape[0])
        print("Difference between Max(std) and Mean(std): ", MaxStdDiffMean[-1])
        print("==============================================================================")

        # SAVE dictionary GPR, samples, difference
        with open(folder+'/Dictionary.pkl','wb') as f:
            pickle.dump({"GPR":gpr, "samples_parameters":samples, "difference_std": MaxStdDiffMean},f)
    else:
        Par = np.zeros(Npar+1, dtype='d')
        while (Par[-1]==0.0):
            comm.Recv(Par, source=0,tag=rank)
            if (Par[-1] == 1.0): break
            OUT = Model(Par[:-1])
            comm.Send(OUT, dest=0,tag=rank+size)
