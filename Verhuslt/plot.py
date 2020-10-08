import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def PlotPosterior(Parameter, color): 
  sns.set()
  sns.set_style('white')
  fig, ax = plt.subplots(1,Parameter.shape[1])
  for i in range(0, Parameter.shape[1]):  
    value = sns.distplot(Parameter[:,i],color=color,ax=ax[i]).get_lines()[0].get_data()
    maxPar = value[0][np.argmax(value[1])]
    ax[i].set_title("MAP = %2.4f" % (maxPar), fontsize=18)
    if (i==0): ax[0].set_xlabel(r'$r$',fontsize=18)
    if (i==1): ax[1].set_xlabel(r'$\kappa$',fontsize=18)
    ax[0].set_ylabel('Density',fontsize=18)
  plt.subplots_adjust(left=0.13,right=0.95,bottom=0.15,top=0.67,wspace=0.38)
  
def Read_Plot(file,color): # reading and plotting paramaters calibrated
  input = np.loadtxt(file, dtype='f', delimiter=' ')
  par1 = np.array(input[:,0])
  par2 = np.array(input[:,1])
  MatrixPar = np.column_stack((par1,par2))
  PlotPosterior(MatrixPar, color)
  
Read_Plot("CalibSMC.dat",color="gray")
Read_Plot("CalibMCMC.dat",color="gray")
Read_Plot("CalibMCMC_AGPR.dat",color="blue")
Read_Plot("CalibSMC_AGPR.dat",color="blue")
plt.show()