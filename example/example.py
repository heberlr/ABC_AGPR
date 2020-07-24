import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
sys.path.append('../src')
from ABC import ABC_SMC,ABC_MCMC
from AGPR import AdapGP

def Reta(Par):
  x = np.arange(0,1,0.01)
  return Par[0]*x + Par[1]
  
def RetaAGPR(Par):
  ParLocal = np.reshape(Par, (1, Par.size))
  with open('model.pkl', 'rb') as f:
    gpr = pickle.load(f)
  OUT = gpr.predict(ParLocal, return_std=True) 
  return OUT[0].ravel()

def CalibReta():
  # Observational data
  x = np.arange(0,1,0.01)
  fx = 1.0*x + 0.5 # The exact parameter is theta = ( 1.0 , 0.5 ) 
  # Boundary of parameter space
  UpperLimit = np.array([1.5,1.0])
  LowLimit = np.array([0.5,0.0])
  # Tolerance vector
  epsilon = np.array([1.0,0.5,0.1])
  # Calling the method SMC
  outSMC = ABC_SMC(Reta, fx, LowLimit, UpperLimit,'CalibSMC.dat',epsilon, 100)
  PlotPosterior(outSMC,([1,2]))
  # Calling the method MCMC
  outMCMC = ABC_MCMC(Reta, fx, LowLimit, UpperLimit,'CalibMCMC.dat',epsilon[-1], 100)
  PlotPosterior(outMCMC,([1,2]))
  # Calling the method SMC-AGPR
  outSMC = ABC_SMC(RetaAGPR, fx, LowLimit, UpperLimit,'CalibSMC.dat',epsilon, 100)
  PlotPosterior(outSMC,([1,2]))
  # Calling the method MCMC-AGPR
  outMCMC = ABC_MCMC(RetaAGPR, fx, LowLimit, UpperLimit,'CalibMCMC.dat',epsilon[-1], 100)
  PlotPosterior(outMCMC,([1,2]))
  
  # Testing prediction of GPR
  # x = np.arange(0,1,0.01) # Same value of the Reta() and RetaGPR()
  # theta = np.array([0.5,1.0]) # Parameter vector, you can change the values and compare them.
  # plt.plot(x, Reta(theta),c='gray',lw=2.5)
  # plt.plot(x, RetaGPR(theta),ls=':',c='red',lw=2.5)
  # plt.xlabel('x')
  # plt.xlabel('f(x)')
  # plt.show()
  
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

def Model(Par):
  run = "./MODEL/replica.sh "+str(Par[0])+" "+str(Par[1])
  #print(run)
  args = shlex.split(run)
  p = subprocess.Popen(args)
  if p.wait() != 0:
    print("There was an error")

  input = np.loadtxt("./MODEL/statistic/confluence-0000.dat", dtype='f', delimiter='\t')
  live = np.array(input[:,1])
  dead = np.array(input[:,3])
  return np.concatenate((live,dead),axis=None)

x = np.arange(0,1,0.01)
# Boundary of parameter space
UpperLimit = np.array([1.5,1.0])
LowLimit = np.array([0.5,0.0])

#Training the Gaussian Process Regression (uncomment to train the GPR)
#AdapGP(Reta,50, LowLimit, UpperLimit, x.shape[0], tol = 1.0)

# Running Reta example
CalibReta()

