import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from verhulst import Model,GenerateData
import scipy.stats as sts
import pickle #library to store the model trained

def CustomPlotPosterior(Parameter, color):
  sns.set()
  sns.set_style('white')
  MAP = np.zeros(Parameter.shape[1])
  fig, ax = plt.subplots(1,Parameter.shape[1])
  for i in range(0, Parameter.shape[1]):
    value = sns.distplot(Parameter[:,i],color=color,ax=ax[i]).get_lines()[0].get_data()
    MAP[i] = value[0][np.argmax(value[1])]
    ax[i].set_title("MAP = %2.4f" % (MAP[i]), fontsize=18)
    if (i==0): ax[0].set_xlabel(r'$r$',fontsize=18)
    if (i==1): ax[1].set_xlabel(r'$\kappa$',fontsize=18)
    ax[0].set_ylabel('Density',fontsize=18)
  plt.subplots_adjust(left=0.13,right=0.95,bottom=0.15,top=0.67,wspace=0.38)
  return MAP

def Read_Plot(file,color): # reading and plotting paramaters calibrated
  input = np.loadtxt(file, dtype='float64', delimiter=' ')
  par1 = np.array(input[:,0])
  par2 = np.array(input[:,1])
  MatrixPar = np.column_stack((par1,par2))
  return CustomPlotPosterior(MatrixPar, color)


def Plot_MAP_response():
  MAP1 = Read_Plot("Calibration/CalibMCMC.dat",color="gray")
  MAP2 = Read_Plot("Calibration/CalibMCMC_AGPR.dat",color="blue")
  MAP3 = Read_Plot("Calibration/CalibSMC.dat",color="gray")
  MAP4 = Read_Plot("Calibration/CalibSMC_AGPR.dat",color="blue")
  # plt.show()

  QOI1 = Model(np.array(MAP1)) # MCMC
  QOI2 = Model(np.array(MAP2)) # MCMC_AGPR
  QOI3 = Model(np.array(MAP3)) # SMC
  QOI4 = Model(np.array(MAP4)) # SMC_AGPR

  data = GenerateData()
  time = np.linspace(0,288,13) #0-12 days

  plt.figure(figsize=(12,8))
  plt.errorbar(time,np.insert(data,0,1000.0),yerr=200.0,fmt='o',color='black',label='Obs. data')
  plt.plot(time,np.insert(QOI1,0,1000.0),color='gray',label='MCMC',ls='-.')
  plt.plot(time,np.insert(QOI2,0,1000.0),color='blue',label='MCMC_AGPR',ls='-.')
  plt.plot(time,np.insert(QOI3,0,1000.0),color='gray',label='SMC',ls=':')
  plt.plot(time,np.insert(QOI4,0,1000.0),color='blue',label='SMC_AGPR',ls=':')
  plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),ncol=6,fontsize=14)
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
  plt.xlabel('Time (hours)', fontsize=18)
  plt.ylabel('Number of cells', fontsize=18)
  plt.subplots_adjust(left=0.14,right=0.90,bottom=0.14,top=0.80)
  plt.show()


def Plot_boxplot(): # plotting quartiles of methods
  def color_box(bp, color):
    # Define the elements to color. You can also add medians, fliers and means
    elements = ['boxes','caps','whiskers','medians']

    # Iterate over each of the elements changing the color
    for elem in elements:
        #[plt.setp(bp[elem][idx], color=color_temp) for idx in range(len(bp[elem]))]
        for idx in range(len(bp[elem])):
          if (idx % 2 == 0 and elem == 'boxes'): color_temp = "gray"
          if (idx % 2 == 1 and elem == 'boxes'): color_temp = "blue"
          if (idx % 4 < 2 and elem == 'caps'): color_temp = "gray"
          if (idx % 4 > 1 and elem == 'caps'): color_temp = "blue"
          if (idx % 4 < 2 and elem == 'whiskers'): color_temp = "gray"
          if (idx % 4 > 1 and elem == 'whiskers'): color_temp = "blue"
          if (elem == 'medians'):
            color_temp = color
            linewidth_temp = 2
          else:
            linewidth_temp = 3
          plt.setp(bp[elem][idx], color=color_temp, linewidth = linewidth_temp)
    return

  boxprops = dict(linestyle='--', linewidth=3, color='darkgoldenrod')
  flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                  linestyle='none')
  medianprops = dict(linestyle='-.', linewidth=2.5, color='firebrick')
  meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick')
  meanlineprops = dict(linestyle='--', linewidth=2.5, color='purple')

  input1 = np.loadtxt("Calibration/CalibMCMC.dat", dtype='float64', delimiter=' ')
  input2 = np.loadtxt("Calibration/CalibSMC.dat", dtype='float64', delimiter=' ')
  input3 = np.loadtxt("Calibration/CalibMCMC_AGPR.dat", dtype='float64', delimiter=' ')
  input4 = np.loadtxt("Calibration/CalibSMC_AGPR.dat", dtype='float64', delimiter=' ')

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
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar1, showfliers=False, labels=["MCMC","MCMC_AGPR","SMC","SMC_AGPR"])
  color_box(bp, 'red')
  plt.ylabel(r'$r$', fontsize=18)
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar2, showfliers=False, labels=["MCMC","MCMC_AGPR","SMC","SMC_AGPR"])
  color_box(bp, 'red')
  plt.ylabel(r'$\kappa$', fontsize=18)
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.11,top=0.88)
  plt.show()

def Kullback_Leibler_divergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def smoothed_hist_kl_distance(fileA, fileB, par, nbins=100, sigma=0.0001):
    inputA = np.loadtxt(fileA, dtype='float64', delimiter=' ')
    inputB = np.loadtxt(fileB, dtype='float64', delimiter=' ')
    parA = np.array(inputA[:,par])
    parB = np.array(inputB[:,par])
    # ahist, bhist = (np.histogram(parA, bins=nbins)[0],np.histogram(parB, bins=nbins)[0])
    # asmooth, bsmooth = (gaussian_filter(ahist, sigma), gaussian_filter(bhist, sigma))
    min = np.min([parA,parB])
    max = np.max([parA,parB])
    DensityA, barA = np.histogram(parA,bins=nbins,range=(min,max),density=True)
    DensityB, barB = np.histogram(parB,bins=nbins,range=(min,max),density=True)
    kdeA = sts.gaussian_kde(parA)
    kdeB = sts.gaussian_kde(parB)
    x = np.linspace(min, max,nbins)
    NormDensityA = DensityA/sum(DensityA)
    NormDensityB = DensityB/sum(DensityB)
    norm_kdeA = kdeA.pdf(x)/sum(kdeA.pdf(x))
    norm_kdeB = kdeB.pdf(x)/sum(kdeB.pdf(x))
    # plt.bar(barA[:-1], NormDensityA, width=np.diff(barA), ec='k', align='edge', label='histogram_A',alpha=0.5)
    # plt.bar(barB[:-1], NormDensityB, width=np.diff(barB), ec='k', align='edge', label='histogram_B',alpha=0.5)
    # plt.plot(x, norm_kdeA, c='C0', lw=4, label='KDE_A')
    # plt.plot(x, norm_kdeB, c='C1', lw=4, label='KDE_B')
    # plt.show()
    return Kullback_Leibler_divergence(np.array(norm_kdeA), np.array(norm_kdeB))



def PlotAGPR_Convergence():
    with open('AGPR/Dictionary.pkl', 'rb') as f:
      dict = pickle.load(f)
    MatrixPar = dict['samples_parameters']
    DiffSTD = dict['difference_std']
    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(12,4))
    iteration = np.arange(0,DiffSTD.shape[0],1)
    ax1.plot(iteration,DiffSTD, marker='.', color='r', markersize=5)
    ax1.set_title('Difference between max(std) and mean(std)')
    ax2.scatter(MatrixPar[:,0],MatrixPar[:,1],s=0.75)
    ax2.set_title('Samples')
    plt.show()


if __name__ == '__main__':
    # Calibration ABC
    Plot_MAP_response()

    # Kullback-Leibler divergence
    print("K-L (SMC,SMC_AGPR) for r: %e" % smoothed_hist_kl_distance("Calibration/CalibSMC.dat", "Calibration/CalibSMC_AGPR.dat", par=0))
    print("K-L (MCMC,MCMC_AGPR) for r: %e" % smoothed_hist_kl_distance("Calibration/CalibMCMC.dat", "Calibration/CalibMCMC_AGPR.dat", par=0))
    print("K-L (SMC,SMC_AGPR) for K: %e" % smoothed_hist_kl_distance("Calibration/CalibSMC.dat", "Calibration/CalibSMC_AGPR.dat", par=1))
    print("K-L (MCMC,MCMC_AGPR) for K: %e" % smoothed_hist_kl_distance("Calibration/CalibMCMC.dat", "Calibration/CalibMCMC_AGPR.dat", par=1))

    # Percentile of distributions
    Plot_boxplot()

    # Plotting convergence of AGPR
    PlotAGPR_Convergence()
