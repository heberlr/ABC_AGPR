import numpy as np
import pickle #library to store the model trained

# Changing repository
import sys
sys.path.append('../src')
# Loading function to ABC calibration
from ABC import ABC_SMC,ABC_MCMC
from AGPR import AdapGP
from GPR import GPR

def GenerateData():
  # Fixed seed - parameters: growth rate = 0.05 and carrying capacity = 10000.0
  rng = np.random.RandomState(1234)
  Par = np.array([0.05,10000.0])
  t = np.linspace(24,288,12) #1-12 days
  OutputModel = np.zeros(t.size)
  InitialCond = 1000.0
  std = InitialCond*0.2
  OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1]) + rng.normal(0, std, t.size)
  return OutputModel

def Model(Par):
  t = np.linspace(24,288,12)
  OutputModel = np.zeros(t.size)
  InitialCond = 1000.0
  OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1])
  return OutputModel

def ModelAGPR(Par):
  ParLocal = np.reshape(Par, (1, Par.size))
  with open('AGPR/Dictionary.pkl', 'rb') as f:
    dict = pickle.load(f)
  OUT = dict['GPR'].predict(ParLocal, return_std=True)
  return OUT[0].ravel()

if __name__ == '__main__':
  # Observational data
  Data = GenerateData()

  # Boundary of parameter space (prior distribution always uniform)
  UpperLimit = np.array([0.1,1.5e4])
  LowLimit = np.array([0.0,0.5e4])
  # Tolerance vector
  epsilon = np.array([1.0e4,0.5e4,0.2e4])
  # Number of accepted parameters that it will generate a posterior distribution
  NumAccept = 1000

  # Calling the method SMC
  outSMC = ABC_SMC(Model, Data, LowLimit, UpperLimit,'Calibration/CalibSMC.dat',epsilon, NumAccept)
  # Calling the method MCMC
  outMCMC = ABC_MCMC(Model, Data, LowLimit, UpperLimit,'Calibration/CalibMCMC.dat',epsilon[-1], NumAccept)

  #Training the Gaussian Process Regression (uncomment to train the GPR)
  NpartionsLHD = 10 # Number of partitions of the Latin hypercube design
  AdapGP(Model,NpartionsLHD, LowLimit, UpperLimit, NumQOI = Data.size, folder='AGPR', tol = 1e-3, DataMesh=True)

  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(ModelAGPR, Data, LowLimit, UpperLimit,'Calibration/CalibSMC_AGPR.dat',epsilon, NumAccept)

  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(ModelAGPR, Data, LowLimit, UpperLimit,'Calibration/CalibMCMC_AGPR.dat',epsilon[-1], NumAccept)
