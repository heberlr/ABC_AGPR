import matplotlib.pyplot as plt
import numpy as np
import pickle

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
  
def differenceSTD_GPR():
    Npar = 2
    NsubDiv = 20
    UpperLimit = np.array([0.5,0.5])
    LowLimit = np.array([0.0,0.0])
    #Mesh for the parametric space
    x_0 = np.zeros(Npar)
    x_f = np.ones(Npar)

    ParMesh = CreateMesh(NsubDiv+1,x_0,x_f) # Parametric space unity
    PosNode = AssPosMesh(ParMesh) # Reference of the nodes associate to elemets
    
    #Transform unit hypercube in real values
    for j in range(0, Npar):
        ParMesh[:,j] = LowLimit[j] + ParMesh[:,j]*(UpperLimit[j]-LowLimit[j])
       
    MaxStdDiffMean = []
    TimeAnalysis = [100,200,300,400,500,600,700,800,900,1000]     
    for ModelSamples in TimeAnalysis:
        #Load model
        fileModel = 'GPR_'+"%03i"%ModelSamples+'/modelGPR.pkl'
        with open(fileModel, 'rb') as f:
            gpr = pickle.load(f)
          
        # Obtain prectidion for the mesh points
        GPstd = np.zeros((NsubDiv+1)**Npar)
        StdElem = np.zeros((NsubDiv)**Npar)
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
    return MaxStdDiffMean 

def Plot_STD(file): # plotting convergence of AGPR
    StdSamples = [0.018475716562097746, 0.004622721762390119, 0.005042477147187073, 0.00515018848693156, 0.005439574263509712, 0.001784377297374005, 0.0014005182004190988, 0.001414468733444876, 0.001392501343749639, 0.004057643966617663]
    Nsamples = [100,200,300,400,500,600,700,800,900,1000]
    input = np.loadtxt(file, dtype='f', delimiter=' ')
    STD = np.array(input)
    it = np.linspace(24,STD.shape[0]+24,STD.shape[0])
    plt.plot(it,STD,"o",color='blue')
    plt.plot(Nsamples,StdSamples,"o",color='red')
    plt.xlabel("Number of samples", fontsize=14)
    plt.ylabel("$log(d_\sigma)$", fontsize=14)
    plt.yscale('log')
    # plt.ylabel("$d_\sigma$ (logarithmic scale)", fontsize=12)
    plt.show()
  
Plot_STD("../Output_AGPR/diffenceSTD.txt")