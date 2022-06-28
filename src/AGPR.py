import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import pickle
from scipy.stats import qmc

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

def ElemSampling(Ind, ParMesh, PosNode, rng):
  N_i = ParMesh[PosNode[Ind,0],:] # lower node which generate the element
  N_f = ParMesh[PosNode[Ind,-1],:] # upper node which generate the element
  Delta = N_f-N_i # vector of differences
  return N_i + (Delta*rng.uniform(0,1,size=Delta.shape[0]))

def AdapGP(Model, NpartionsLHD, LowLimit, UpperLimit, NumQOI, folder, max_iterations=10000, tol=10**-3, NSampCov=40, DataMesh=False, SeedRandomState = 1234):
  rng = np.random.RandomState(SeedRandomState) # random number generator
  Npar = UpperLimit.shape[0]
  lhs = qmc.LatinHypercube(d=Npar,seed=rng)
  samples = lhs.random(n=NpartionsLHD) # Latin Hypercube Sampling of Npar variables, and NpartionsLHD samples each.
  vertices = CreateVert(Npar) # Vertices of the hypercube

  #Concatenate samples and nodes of parametric hypercube
  samples = np.concatenate((samples, vertices), axis=0)

  #Mesh for the parametric space
  x_0 = np.zeros(Npar)
  x_f = np.ones(Npar)

  ParMesh = CreateMesh(NpartionsLHD+1,x_0,x_f) # Parametric space unity (shape = (number of nodes, number of parameters))
  PosNode = AssPosMesh(ParMesh) # Reference of the nodes associate to elemets (shape = (number of elements, number of connections by node))
  dict_MeshValues = {}

  #Transform unit hypercube in real values
  for j in range(0, Npar):
    samples[:,j] = LowLimit[j] + samples[:,j]*(UpperLimit[j]-LowLimit[j])
    ParMesh[:,j] = LowLimit[j] + ParMesh[:,j]*(UpperLimit[j]-LowLimit[j])

  #Output of the model
  OutputModel = np.zeros( (samples.shape[0],NumQOI) )
  for i in range(0, samples.shape[0]):
    OutputModel[i,:] = Model(samples[i,:])

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
  if (DataMesh):
    dict_MeshValues['Mesh'] = ParMesh
    dict_MeshValues['Nodes_Element'] = PosNode
    dict_MeshValues[0] = {'GPstdNodes': GPstd.copy(), 'StdElements': StdElem.copy()}

  #Calculating difference between StdMax amd StdMean
  MeanStd = np.mean(StdElem)
  InitialSamples = samples
  MaxStdDiffMean = []
  MaxStdDiffMean.append(StdElem.max() - MeanStd)
  print("Initial step...")
  # Main loop
  for k in range(0,max_iterations):
    Ind = np.argmax(StdElem)
    value = ElemSampling(Ind,ParMesh,PosNode,rng)
    AddSolution = np.array([Model(value)])
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
    if (DataMesh):
        dict_MeshValues[k+1] = {'GPstdNodes': GPstd.copy(), 'StdElements': StdElem.copy()}

    #Calculating difference between StdMax amd StdMean
    MeanStd = np.mean(StdElem)
    MaxStdDiffMean.append(StdElem.max() - MeanStd)

    print("Number of new samples: "+str(k+1)+" Crit_Diff: "+str(MaxStdDiffMean[-1]))
    #Stop criterion
    if (k >= NSampCov) and (CritConv(MaxStdDiffMean[-NSampCov:]) < tol):
      break

  #Finishing the metamodel
  print("==============================================================================")
  print("Total of iterations: ", k+1)
  print("Total of samples: ", samples.shape[0])
  print("Difference between Max(std) and Mean(std): ", MaxStdDiffMean[-1])
  print("==============================================================================")

  # SAVE dictionary GPR, samples, difference
  with open(folder+'/Dictionary.pkl','wb') as f:
      pickle.dump({"GPR":gpr, "samples_parameters":samples, "difference_std": MaxStdDiffMean},f)
  if (DataMesh):
      with open(folder+'/MeshValues.pkl','wb') as f:
          pickle.dump(dict_MeshValues,f)
