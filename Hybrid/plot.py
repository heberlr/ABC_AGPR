import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def PlotPosterior(Parameter, color): 
  fig = plt.figure()
  sns.set_color_codes("deep")
  sns.set_style('darkgrid')
  fig.subplots_adjust(hspace=0.4, wspace=0.4)
  for i in range(1, Parameter.shape[1]+1):
    ax = fig.add_subplot(1, Parameter.shape[1], i)
    #value = sns.distplot(Parameter[:,i-1],color='blue',bins=10).get_lines()[0].get_data()
    value = sns.distplot(Parameter[:,i-1],color=color).get_lines()[0].get_data()
    maxPar = value[0][np.argmax(value[1])]
    plt.xlabel(r"Parameter $\bar{\alpha}_{P}$ - MAP = %2.4f" % (maxPar), fontsize=12)
  
def Read_Plot(file,color):
  input = np.loadtxt(file, dtype='f', delimiter=' ')
  par1 = np.array(input[:,0])
  par2 = np.array(input[:,1])
  MatrixPar = np.column_stack((par1,par2))
  PlotPosterior(MatrixPar, color)
  
def Plot_STD(file):
  input = np.loadtxt(file, dtype='f', delimiter=' ')
  STD = np.array(input)
  it = np.linspace(0,STD.shape[0],STD.shape[0])
  plt.plot(it,STD,"o")
  plt.xlabel("Number of iterations", fontsize=14)
  plt.ylabel("$d_\sigma$", fontsize=14)
  # plt.yscale('log')
  # plt.ylabel("$d_\sigma$ (logarithmic scale)", fontsize=12)
  plt.show()

def Plot_boxplot():
  def color_box(bp, color):

    # Define the elements to color. You can also add medians, fliers and means
    elements = ['boxes','caps','whiskers','medians']

    # Iterate over each of the elements changing the color
    for elem in elements:
        #[plt.setp(bp[elem][idx], color=color_temp) for idx in range(len(bp[elem]))]
        for idx in range(len(bp[elem])):
          if (idx % 2 == 0 and elem == 'boxes'): color_temp = "blue"
          if (idx % 2 == 1 and elem == 'boxes'): color_temp = "black"
          if (idx % 4 < 2 and elem == 'caps'): color_temp = "blue"
          if (idx % 4 > 1 and elem == 'caps'): color_temp = "black"
          if (idx % 4 < 2 and elem == 'whiskers'): color_temp = "blue"
          if (idx % 4 > 1 and elem == 'whiskers'): color_temp = "black"
          if (elem == 'medians'): color_temp = color
          plt.setp(bp[elem][idx], color=color_temp)
    return
    
  input1 = np.loadtxt("CalibMCMC.dat", dtype='f', delimiter=' ')
  input2 = np.loadtxt("CalibSMC.dat", dtype='f', delimiter=' ')
  input3 = np.loadtxt("CalibMCMC_AGPR.dat", dtype='f', delimiter=' ')
  input4 = np.loadtxt("CalibSMC_AGPR.dat", dtype='f', delimiter=' ')

  Inp1Par1 = np.array(input1[:,0])
  Inp1Par2 = np.array(input1[:,1])
  Inp2Par1 = np.array(input2[:,0])
  Inp2Par2 = np.array(input2[:,1])
  Inp3Par1 = np.array(input3[:,0])
  Inp3Par2 = np.array(input3[:,1])
  Inp4Par1 = np.array(input4[:,0])
  Inp4Par2 = np.array(input4[:,1])
  MatrixPar1 = np.column_stack((Inp1Par1,Inp3Par1,Inp2Par1,Inp4Par1))
  MatrixPar2 = np.column_stack((Inp1Par2,Inp3Par2,Inp2Par2,Inp4Par2))
  plt.figure()
  bp = plt.boxplot(MatrixPar1, showfliers=False, labels=["MCMC","MCMC_AGPR","SMC","SMC_AGPR"])
  color_box(bp, 'red')
  plt.ylabel(r'Parameter $\bar{\alpha}_{P}$', fontsize=12)
  plt.figure()
  bp = plt.boxplot(MatrixPar2, showfliers=False, labels=["MCMC","MCMC_AGPR","SMC","SMC_AGPR"])
  color_box(bp, 'red')
  plt.ylabel(r'Parameter $\bar{\alpha}_{D}$', fontsize=12)
  plt.show()

# Read_Plot("CalibSMC.dat",color="gray")
# Read_Plot("CalibMCMC.dat",color="gray")
# Read_Plot("CalibMCMC_AGPR.dat",color="blue")
# Read_Plot("CalibSMC_AGPR.dat",color="blue")
# plt.show()  
# Plot_STD("diffenceSTD.txt")
# Plot_boxplot()