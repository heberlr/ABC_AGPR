import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sts
import pickle #library to store the model trained

def ReadParameters(File):
    input = np.loadtxt(File, dtype='float64', delimiter=' ')
    par1 = np.array(input[:,0])
    par2 = np.array(input[:,1])
    return np.column_stack((par1,par2))

def Get_MAP(Parameter):
    MAP = np.zeros(Parameter.shape[1])
    for i in range(0, Parameter.shape[1]):
        plt.clf()
        value = sns.kdeplot(Parameter[:,i]).get_lines()[0].get_data()
        MAP[i] = value[0][np.argmax(value[1])]
        plt.title("MAP = %2.4f" % (MAP[i]), fontsize=18)
        #plt.axvline(MAP[i],color='Red')
        #plt.show()
        plt.clf()
        plt.close()
    return MAP

def PlotPosterior(file, color): # plot histograms with fitting of the density curves
  Parameter = ReadParameters(file)
  MAP = Get_MAP(Parameter)
  # Plotting
  sns.set()
  sns.set_style('ticks')
  fig, ax = plt.subplots(1,Parameter.shape[1])
  for i in range(0, Parameter.shape[1]):
    sns.histplot(Parameter[:,i],stat='probability',color=color, kde=True, ax=ax[i]) # bins=nbins,
    ax[i].axvline(MAP[i],color='Red')
    ax[i].set_title("MAP = %2.4f" % (MAP[i]), fontsize=18)
    if (i==0): ax[0].set_xlabel(r'$\bar{\alpha}_{P}$',fontsize=18)
    if (i==1): ax[1].set_xlabel(r'$\bar{\alpha}_{D}$',fontsize=18)
  ax[0].set(ylim = (0,0.13),xticks=(np.linspace(0,0.5,3)))
  ax[1].set(ylabel=None,ylim = (0,0.12), xticks=(np.linspace(0.1,0.25,4)))
  plt.subplots_adjust(left=0.13,right=0.95,bottom=0.15,top=0.67,wspace=0.38)
  FileOut = file.split("/")[-1].split(".")[0]+".jpg"
  fig.savefig(FileOut, dpi=120,bbox_inches='tight')

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
    fig.savefig("AGPR_iterations.jpg", dpi=120,bbox_inches='tight')

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

  input1 = np.loadtxt("./Calibration/CalibMCMC.dat", dtype='float64', delimiter=' ')
  input2 = np.loadtxt("./Calibration/CalibSMC.dat", dtype='float64', delimiter=' ')
  input3 = np.loadtxt("./Calibration/CalibMCMC_AGPR.dat", dtype='float64', delimiter=' ')
  input4 = np.loadtxt("./Calibration/CalibSMC_AGPR.dat", dtype='float64', delimiter=' ')

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
  plt.ylabel(r'$\bar{\alpha}_{P}$', fontsize=18)
  plt.savefig("BoxPlot_alpha_P.jpg", dpi=120,bbox_inches='tight')
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar2, showfliers=False, labels=["MCMC","MCMC_AGPR","SMC","SMC_AGPR"])
  color_box(bp, 'red')
  plt.ylabel(r'$\bar{\alpha}_{D}$', fontsize=18)
  plt.savefig("BoxPlot_alpha_D.jpg", dpi=120,bbox_inches='tight')
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.11,top=0.88)

def Plot_MAP_response():
  input1 = np.loadtxt("./Output_SingleRun/OutModel_ABC_MCMC.dat", dtype='float64', delimiter='\t')
  input2 = np.loadtxt("./Output_SingleRun/OutModel_ABC_MCMC_AGPR.dat", dtype='float64', delimiter='\t')
  input3 = np.loadtxt("./Output_SingleRun/OutModel_ABC_SMC.dat", dtype='float64', delimiter='\t')
  input4 = np.loadtxt("./Output_SingleRun/OutModel_ABC_SMC_AGPR.dat", dtype='float64', delimiter='\t')
  QOI1 = np.array(input1)
  QOI2 = np.array(input2)
  QOI3 = np.array(input3)
  QOI4 = np.array(input4)
  Mean1 = np.mean(QOI1,axis=0)
  Mean2 = np.mean(QOI2,axis=0)
  Mean3 = np.mean(QOI3,axis=0)
  Mean4 = np.mean(QOI4,axis=0)
  data = np.array([0.56829, 0.5681037, 0.5686626, 0.5691672, 0.6176569, 0.6802602, 0.7250002, 0.7600671, 0.7605562, 0.7590967, 0.7202258, 0.6275473, 0.538145, 0.4514444, 0.3496518, 0.2844012, 0.1959615, 0.1393592, 0.1014277, 0.06882953, 0.04892439, 0.0309368, 0.02281637, 0.01860866, 0.01670665, 0.01145865, 0.00999915, 0.01005349, 0.00989046, 0.00784871, 0.0077245, 0.00554301, 0.00337704, 0.01571294, 0.0155732, 0.0155111, 0.01561202, 0.0151928, 0.01500648, 0.01520833, 0.0155111, 0.01558873, 0.01829036, 0.07083246, 0.1857993, 0.2866292, 0.3938561, 0.4960524, 0.5698349, 0.6533526, 0.7113057, 0.7336408, 0.7633044, 0.7718518, 0.784506, 0.7873551, 0.7897307, 0.7894745, 0.7897229, 0.7902896, 0.7875647, 0.786703, 0.785818, 0.784638, 0.7825574, 0.7837685 ])

  time = np.linspace(0,96,33)
  plt.figure(figsize=(12,8))
  plt.errorbar(time,data[0:33],yerr=0.03,fmt='o',color='green',label='Live')
  plt.plot(time,Mean1[0:33],color='gray',label='MCMC',ls='-.')
  plt.plot(time,Mean2[0:33],color='blue',label='MCMC_AGPR',ls='-.')
  plt.plot(time,Mean3[0:33],color='gray',label='SMC',ls=':')
  plt.plot(time,Mean4[0:33],color='blue',label='SMC_AGPR',ls=':')
  plt.errorbar(time,data[33:],yerr=0.03,fmt='o',color='red',label='Dead')
  plt.plot(time,Mean1[33:],color='gray',ls='-.')
  plt.plot(time,Mean2[33:],color='blue',ls='-.')
  plt.plot(time,Mean3[33:],color='gray',ls=':')
  plt.plot(time,Mean4[33:],color='blue',ls=':')
  plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),ncol=6,fontsize=14)
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
  plt.xlabel('Time (hours)', fontsize=18)
  plt.ylabel('Confluence', fontsize=18)
  plt.subplots_adjust(left=0.14,right=0.90,bottom=0.14,top=0.80)
  plt.savefig("MAP_response.jpg", dpi=120,bbox_inches='tight')

def Kullback_Leibler_divergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def smoothed_hist_kl_distance(fileA, fileB, par, nbins=20, sigma=0.0001):
    inputA = np.loadtxt(fileA, dtype='float64', delimiter=' ')
    inputB = np.loadtxt(fileB, dtype='float64', delimiter=' ')
    parA = np.array(inputA[:,par])
    parB = np.array(inputB[:,par])
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
    return Kullback_Leibler_divergence(np.array(norm_kdeA), np.array(norm_kdeB))

if __name__ == '__main__':
    # # Calibration ABC
    PlotPosterior("./Calibration/CalibSMC.dat",color="gray")
    PlotPosterior("./Calibration/CalibMCMC.dat",color="gray")
    # Calibration ABC with AGPR
    PlotPosterior("./Calibration/CalibMCMC_AGPR.dat",color="blue")
    PlotPosterior("./Calibration/CalibSMC_AGPR.dat",color="blue")

    # Kullback-Leibler divergence
    print("K-L (SMC,SMC_AGPR) for alpha_P: %e" % smoothed_hist_kl_distance("./Calibration/CalibSMC.dat", "./Calibration/CalibSMC_AGPR.dat", par=0))
    print("K-L (MCMC,MCMC_AGPR) for alpha_P: %e" % smoothed_hist_kl_distance("./Calibration/CalibMCMC.dat", "./Calibration/CalibMCMC_AGPR.dat", par=0))
    print("K-L (SMC,SMC_AGPR) for alpha_D: %e" % smoothed_hist_kl_distance("./Calibration/CalibSMC.dat", "./Calibration/CalibSMC_AGPR.dat", par=1))
    print("K-L (MCMC,MCMC_AGPR) for alpha_D: %e" % smoothed_hist_kl_distance("./Calibration/CalibMCMC.dat", "./Calibration/CalibMCMC_AGPR.dat", par=1))

    # Percentile of distributions
    Plot_boxplot()

    # Simulation of MAPs from each posterior distribution
    Plot_MAP_response()

    # Plotting convergence of AGPR
    PlotAGPR_Convergence()
