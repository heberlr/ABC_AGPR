import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
sys.path.append('../src')
from ABC import ABC_SMC,ABC_MCMC
from AGPR import AdapGP

def GenerateData():
  # Fixed seed
  np.random.seed(1234)
  Par = np.array([0.05,10000.0]) 
  t = np.linspace(24,288,12)
  OutputModel = np.zeros(t.size)
  InitialCond = 1000.0
  variance = InitialCond*0.5
  OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1]) + np.random.normal(0, variance, t.size)
  plt.plot(np.concatenate((0, t), axis=None),np.concatenate((InitialCond, OutputModel), axis=None),"o")
  plt.show()
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
  t = np.linspace(24,288,12)
  Data = GenerateData()
  # Boundary of parameter space
  UpperLimit = np.array([0.1x,1e4])
  LowLimit = np.array([0.0,0.8e4])
  # Tolerance vector
  epsilon = np.array([4.0e4,1.0e4,3.5e3])
  # Calling the method SMC
  outSMC = ABC_SMC(Model, Data, LowLimit, UpperLimit,'CalibSMC.dat',epsilon, 100)
  PlotPosterior(outSMC,([1,2]))
  # Calling the method MCMC
  outMCMC = ABC_MCMC(Model, Data, LowLimit, UpperLimit,'CalibMCMC.dat',epsilon[-1], 100)
  PlotPosterior(outMCMC,([1,2]))
  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(ModelAGPR, Data, LowLimit, UpperLimit,'AGPRCalibSMC.dat',epsilon, 100)
  PlotPosterior(outSMC,([1,2]))
  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(ModelAGPR, Data, LowLimit, UpperLimit,'AGPRCalibMCMC.dat',epsilon[-1], 100)
  PlotPosterior(outMCMC,([1,2]))
  
def PlotPosterior(Parameter, PlotMatrix): 
  fig = plt.figure()
  sns.set_color_codes("deep")
  sns.set_style('darkgrid')
  fig.subplots_adjust(hspace=0.4, wspace=0.4)
  for i in range(1, Parameter.shape[1]+1):
    ax = fig.add_subplot(PlotMatrix[0], PlotMatrix[1], i)
    value = sns.distplot(Parameter[:,i-1],color='blue').get_lines()[0].get_data()
    maxPar = value[0][np.argmax(value[1])]
    plt.xlabel("Parameter %d - MAP = %2.4f" % (i,maxPar), fontsize=12)
  plt.show()

print(GenerateData())
# Boundary of parameter space
#UpperLimit = np.array([0.3,3.7e4])
#LowLimit = np.array([0.0,3.0e4])

#Training the Gaussian Process Regression (uncomment to train the GPR)
#AdapGP(Model,10, LowLimit, UpperLimit, t.size, tol = 1e-3)

#CalibVerhulst()

# Plotting
# plt.errorbar(t, GenerateData(), yerr=1700.0*0.5,marker=".",ls='none')
# plt.plot(t,Model([0.04,35000.0]))
# plt.show()