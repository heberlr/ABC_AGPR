import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
sys.path.append('../src')
from ABC import ABC_SMC,ABC_MCMC
from AGPR import AdapGP

def GenerateData():
  # Fixed seed - parameters: growth rate = 0.05 and carrying capacity = 10000.0
  np.random.seed(1234)
  Par = np.array([0.05,10000.0]) 
  t = np.linspace(24,288,12)
  OutputModel = np.zeros(t.size)
  InitialCond = 1000.0
  std = InitialCond*0.1
  OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1]) + np.random.normal(0, std, t.size)
  return OutputModel

def Model(Par):  
  t = np.linspace(24,288,12)
  OutputModel = np.zeros(t.size)
  InitialCond = 1000.0
  OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1])
  return OutputModel
  
def ModelAGPR(Par):
  ParLocal = np.reshape(Par, (1, Par.size))
  with open('model.pkl', 'rb') as f:
    gpr = pickle.load(f)
  OUT = gpr.predict(ParLocal, return_std=True) 
  return OUT[0].ravel()

def CalibVerhulst():
  # Observational data
  t = np.linspace(24,288,12) #1-12 days
  Data = GenerateData()
  # Boundary of parameter space
  UpperLimit = np.array([0.1,1.5e4])
  LowLimit = np.array([0.0,0.5e4])
  # Tolerance vector
  epsilon = np.array([1.0e4,0.5e4,0.2e4])
  # Calling the method SMC
  outSMC = ABC_SMC(Model, Data, LowLimit, UpperLimit,'CalibSMC.dat',epsilon, 1000)
  # Calling the method MCMC
  outMCMC = ABC_MCMC(Model, Data, LowLimit, UpperLimit,'CalibMCMC.dat',epsilon[-1], 1000)
  
  #Training the Gaussian Process Regression (uncomment to train the GPR)
  #NsamplesInitial=20
  #AdapGP(Model,NsamplesInitial, LowLimit, UpperLimit, t.size, tol = 1e-3)
  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(ModelAGPR, Data, LowLimit, UpperLimit,'CalibSMC_AGPR.dat',epsilon, 1000)
  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(ModelAGPR, Data, LowLimit, UpperLimit,'CalibMCMC_AGPR.dat',epsilon[-1], 1000)


# CalibVerhulst()