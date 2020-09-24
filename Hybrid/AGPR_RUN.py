import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
sys.path.append('../src')
from ABC import ABC_SMC,ABC_MCMC
from AGPR import AdapGP

def Model_AGPR(Par):
  ParLocal = np.reshape(Par, (1, Par.size))
  with open('model.pkl', 'rb') as f:
    gpr = pickle.load(f)
  OUT = gpr.predict(ParLocal, return_std=True) 
  return OUT[0].ravel()

def Calib():
  # Observational data
  data = np.array([0.56829, 0.5681037, 0.5686626, 0.5691672, 0.6176569, 0.6802602, 0.7250002, 0.7600671, 0.7605562, 0.7590967, 0.7202258, 0.6275473, 0.538145, 0.4514444, 0.3496518, 0.2844012, 0.1959615, 0.1393592, 0.1014277, 0.06882953, 0.04892439, 0.0309368, 0.02281637, 0.01860866, 0.01670665, 0.01145865, 0.00999915, 0.01005349, 0.00989046, 0.00784871, 0.0077245, 0.00554301, 0.00337704, 0.01571294, 0.0155732, 0.0155111, 0.01561202, 0.0151928, 0.01500648, 0.01520833, 0.0155111, 0.01558873, 0.01829036, 0.07083246, 0.1857993, 0.2866292, 0.3938561, 0.4960524, 0.5698349, 0.6533526, 0.7113057, 0.7336408, 0.7633044, 0.7718518, 0.784506, 0.7873551, 0.7897307, 0.7894745, 0.7897229, 0.7902896, 0.7875647, 0.786703, 0.785818, 0.784638, 0.7825574, 0.7837685 ]) # [Live; Dead] -> parameter 0.2 and 0.15 with normal(0;0.03) 
  # Boundary of parameter space
  UP = np.array([0.5,0.5])
  DW = np.array([0.0,0.0])
  # Tolerance vector
  epsilon = np.array([1.0,0.5,0.2])
  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(Model_AGPR, data, DW, UP,'CalibSMC_AGPR.dat',epsilon, NumAccept=1000)
  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(Model_AGPR, data, DW, UP,'CalibMCMC_AGPR.dat',epsilon[-1], NumAccept=1000,var_trasition=0.1)
  
  # Testing prediction of GPR
  # x = np.arange(0,1,0.01) # Same value of the Reta() and RetaGPR()
  # theta = np.array([0.5,1.0]) # Parameter vector, you can change the values and compare them.
  # plt.plot(x, Reta(theta),c='gray',lw=2.5)
  # plt.plot(x, RetaGPR(theta),ls=':',c='red',lw=2.5)
  # plt.xlabel('x')
  # plt.xlabel('f(x)')
  # plt.show()
  

# Boundary of parameter space


# Running Reta example
Calib()

