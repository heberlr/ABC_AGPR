import numpy as np
import pickle #library to store the model trained

# Changing repository
import sys
sys.path.append('../src')
# Loading function to ABC calibration
from ABC_parallel import ABC_SMC,ABC_MCMC
from AGPR_parallel import AdapGP
from GPR_parallel import GPR

from mpi4py import MPI
import subprocess
import os

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def ObservationalData():
  data = np.array([0.56829, 0.5681037, 0.5686626, 0.5691672, 0.6176569, 0.6802602, 0.7250002, 0.7600671, 0.7605562, 0.7590967, 0.7202258, 0.6275473, 0.538145, 0.4514444, 0.3496518, 0.2844012, 0.1959615, 0.1393592, 0.1014277, 0.06882953, 0.04892439, 0.0309368, 0.02281637, 0.01860866, 0.01670665, 0.01145865, 0.00999915, 0.01005349, 0.00989046, 0.00784871, 0.0077245, 0.00554301, 0.00337704, 0.01571294, 0.0155732, 0.0155111, 0.01561202, 0.0151928, 0.01500648, 0.01520833, 0.0155111, 0.01558873, 0.01829036, 0.07083246, 0.1857993, 0.2866292, 0.3938561, 0.4960524, 0.5698349, 0.6533526, 0.7113057, 0.7336408, 0.7633044, 0.7718518, 0.784506, 0.7873551, 0.7897307, 0.7894745, 0.7897229, 0.7902896, 0.7875647, 0.786703, 0.785818, 0.784638, 0.7825574, 0.7837685 ])
  return data

def Calling_modelOpenMP(parameter):
    function_call = ['./build/main.exe', '{}'.format(rank)] #change of .exe to .out
    for ind in range(0,parameter.shape[0]):
        function_call.append('{}'.format(parameter[ind]))
    cache = subprocess.run(function_call,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # print(cache.stdout) # print(cache.stderr) # print(cache.returncode)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)
    OUT = np.fromstring(cache.stdout, dtype='d', sep=' ')
    time = OUT[::3]
    Live = OUT[1::3]
    Dead = OUT[2::3]
    NumQOI = 66 # Number of QOIs
    QOI = np.concatenate((Live, Dead), axis=None)
    if (QOI.shape[0] != NumQOI):
        print("Model output error! incompatibility of QoIs!")
        os._exit(0)
    return QOI

def ModelAGPR(Par):
  ParLocal = np.reshape(Par, (1, Par.size))
  with open('AGPR/Dictionary.pkl', 'rb') as f:
    dict = pickle.load(f)
  OUT = dict['GPR'].predict(ParLocal, return_std=True)
  return OUT[0].ravel()

if __name__ == '__main__':
  # Observational data
  Data = ObservationalData()
  # Boundary of parameter space (prior distribution always uniform)
  UpperLimit = np.array([0.5,0.5])
  LowLimit = np.array([0.0,0.0])
  # Tolerance vector
  epsilon = np.array([1.0,0.5,0.2])
  # Number of accepted parameters that it will generate a posterior distribution
  NumAccept = 1000
  # Number of replicates
  Nrep = (size - 1)

  # Calling the method SMC
  outSMC = ABC_SMC(Calling_modelOpenMP, Data, LowLimit, UpperLimit,'Calibration/CalibSMC.dat', Nrep, epsilon, NumAccept)
  # Calling the method MCMC
  outMCMC = ABC_MCMC(Calling_modelOpenMP, Data, LowLimit, UpperLimit,'Calibration/CalibMCMC.dat', Nrep, epsilon[-1], NumAccept)

  #Training the Gaussian Process Regression (uncomment to train the GPR)
  NpartionsLHD = 20 # Number of partitions of the Latin hypercube design
  AdapGP(Calling_modelOpenMP,NpartionsLHD, LowLimit, UpperLimit, NumQOI = Data.size, folder='AGPR', NumRep = Nrep, tol = 1e-3)

  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(Calling_modelOpenMP, Data, LowLimit, UpperLimit,'Calibration/CalibSMC_AGPR.dat', Nrep, epsilon, NumAccept)

  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(Calling_modelOpenMP, Data, LowLimit, UpperLimit,'Calibration/CalibMCMC_AGPR.dat', Nrep, epsilon[-1], NumAccept)
