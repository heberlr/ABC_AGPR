import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def PlotPosterior(Parameter, color): # plot histograms with fitting of the density curves
  sns.set()
  sns.set_style('white')
  fig, ax = plt.subplots(1,Parameter.shape[1])
  for i in range(0, Parameter.shape[1]):  
    value = sns.distplot(Parameter[:,i],color=color,ax=ax[i]).get_lines()[0].get_data()
    maxPar = value[0][np.argmax(value[1])]
    ax[i].set_title("MAP = %2.4f" % (maxPar), fontsize=18)
    if (i==0): 
        ax[0].set_xlabel(r'$\bar{\alpha}_{P}$',fontsize=18)
        ax[0].set_xlim([0,0.6])
    if (i==1): 
        ax[1].set_xlabel(r'$\bar{\alpha}_{D}$',fontsize=18)
        ax[1].set_xlim([0,0.4])
    ax[0].set_ylabel('Density',fontsize=18)
  plt.subplots_adjust(left=0.13,right=0.95,bottom=0.15,top=0.67,wspace=0.38)
  
def Read_Plot(file,color,name): # reading and plotting paramaters calibrated
  input = np.loadtxt(file, dtype='f', delimiter=' ')
  par1 = np.array(input[:,0])
  par2 = np.array(input[:,1])
  MatrixPar = np.column_stack((par1,par2))
  PlotPosterior(MatrixPar, color)
  plt.savefig(name+'.png')
  plt.savefig(name+'.svg')
  
def Plot_STD(file): # plotting convergence of AGPR
  input = np.loadtxt(file, dtype='f', delimiter=' ')
  STD = np.array(input)
  it = np.linspace(0,STD.shape[0],STD.shape[0])
  plt.plot(it,STD,"o")
  plt.xlabel("Number of iterations", fontsize=14)
  plt.ylabel("$log(d_\sigma)$", fontsize=14)
  plt.yscale('log')
  # plt.ylabel("$d_\sigma$ (logarithmic scale)", fontsize=12)
  plt.show()

def Plot_boxplot(): # plotting quartiles of methods
  def color_box(bp, color):

    # Define the elements to color. You can also add medians, fliers and means
    elements = ['boxes','caps','whiskers','medians']

    # Iterate over each of the elements changing the color
    for elem in elements:
        #[plt.setp(bp[elem][idx], color=color_temp) for idx in range(len(bp[elem]))]
        for idx in range(len(bp[elem])):
          # if (idx % 2 == 0 and elem == 'boxes'): color_temp = "gray"
          # if (idx % 2 == 1 and elem == 'boxes'): color_temp = "blue"
          # if (idx % 4 < 2 and elem == 'caps'): color_temp = "gray"
          # if (idx % 4 > 1 and elem == 'caps'): color_temp = "blue"
          # if (idx % 4 < 2 and elem == 'whiskers'): color_temp = "gray"
          # if (idx % 4 > 1 and elem == 'whiskers'): color_temp = "blue"
          if (elem == 'whiskers' or elem == 'caps' or elem == 'boxes'):color_temp = "blue"
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
    
  input1 = np.loadtxt("GPR_100/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input2 = np.loadtxt("GPR_100/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input3 = np.loadtxt("GPR_200/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input4 = np.loadtxt("GPR_200/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input5 = np.loadtxt("GPR_300/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input6 = np.loadtxt("GPR_300/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input7 = np.loadtxt("GPR_400/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input8 = np.loadtxt("GPR_400/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input9 = np.loadtxt("GPR_500/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input10 = np.loadtxt("GPR_500/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input11 = np.loadtxt("GPR_600/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input12 = np.loadtxt("GPR_600/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input13 = np.loadtxt("GPR_700/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input14 = np.loadtxt("GPR_700/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input15 = np.loadtxt("GPR_800/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input16 = np.loadtxt("GPR_800/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input17 = np.loadtxt("GPR_900/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input18 = np.loadtxt("GPR_900/CalibSMC_GPR.dat", dtype='f', delimiter=' ')
  input19 = np.loadtxt("GPR_1000/CalibMCMC_GPR.dat", dtype='f', delimiter=' ')
  input20 = np.loadtxt("GPR_1000/CalibSMC_GPR.dat", dtype='f', delimiter=' ')

  Inp1Par1 = np.array(input1[:,0])
  Inp1Par2 = np.array(input1[:,1])
  Inp2Par1 = np.array(input2[:,0])
  Inp2Par2 = np.array(input2[:,1])
  Inp3Par1 = np.array(input3[:,0])
  Inp3Par2 = np.array(input3[:,1])
  Inp4Par1 = np.array(input4[:,0])
  Inp4Par2 = np.array(input4[:,1])
  Inp5Par1 = np.array(input5[:,0])
  Inp5Par2 = np.array(input5[:,1])
  Inp6Par1 = np.array(input6[:,0])
  Inp6Par2 = np.array(input6[:,1])
  Inp7Par1 = np.array(input7[:,0])
  Inp7Par2 = np.array(input7[:,1])
  Inp8Par1 = np.array(input8[:,0])
  Inp8Par2 = np.array(input8[:,1])
  Inp9Par1 = np.array(input9[:,0])
  Inp9Par2 = np.array(input9[:,1])
  Inp10Par1 = np.array(input10[:,0])
  Inp10Par2 = np.array(input10[:,1])
  Inp11Par1 = np.array(input11[:,0])
  Inp11Par2 = np.array(input11[:,1])
  Inp12Par1 = np.array(input12[:,0])
  Inp12Par2 = np.array(input12[:,1])
  Inp13Par1 = np.array(input13[:,0])
  Inp13Par2 = np.array(input13[:,1])
  Inp14Par1 = np.array(input14[:,0])
  Inp14Par2 = np.array(input14[:,1])
  Inp15Par1 = np.array(input15[:,0])
  Inp15Par2 = np.array(input15[:,1])
  Inp16Par1 = np.array(input16[:,0])
  Inp16Par2 = np.array(input16[:,1])
  Inp17Par1 = np.array(input17[:,0])
  Inp17Par2 = np.array(input17[:,1])
  Inp18Par1 = np.array(input18[:,0])
  Inp18Par2 = np.array(input18[:,1])
  Inp19Par1 = np.array(input19[:,0])
  Inp19Par2 = np.array(input19[:,1])
  Inp20Par1 = np.array(input20[:,0])
  Inp20Par2 = np.array(input20[:,1])
  MatrixPar1 = np.column_stack((Inp1Par1,Inp3Par1,Inp5Par1,Inp7Par1,Inp9Par1,Inp11Par1,Inp13Par1,Inp15Par1,Inp17Par1,Inp19Par1))
  MatrixPar2 = np.column_stack((Inp1Par2,Inp3Par2,Inp5Par2,Inp7Par2,Inp9Par2,Inp11Par2,Inp13Par2,Inp15Par2,Inp17Par2,Inp19Par2))
  MatrixPar3 = np.column_stack((Inp2Par1,Inp4Par1,Inp6Par1,Inp8Par1,Inp10Par1,Inp12Par1,Inp14Par1,Inp16Par1,Inp18Par1,Inp20Par1))
  MatrixPar4 = np.column_stack((Inp2Par2,Inp4Par2,Inp6Par2,Inp8Par2,Inp10Par2,Inp12Par2,Inp14Par2,Inp16Par2,Inp18Par2,Inp20Par2))
  # MCMC alpha_P
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar1, showfliers=False, labels=["100","200","300","400","500","600","700","800","900","1000"])
  color_box(bp, 'red')
  plt.ylabel(r'$\bar{\alpha}_{P}$', fontsize=18)
  plt.xlabel('Number of samples', fontsize=18)
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.12,top=0.88)
  # MCMC alpha_D
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar2, showfliers=False, labels=["100","200","300","400","500","600","700","800","900","1000"])
  color_box(bp, 'red')
  plt.ylabel(r'$\bar{\alpha}_{D}$', fontsize=18)
  plt.xlabel('Number of samples', fontsize=18)
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.12,top=0.88)
  # SMC alpha_P
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar3, showfliers=False, labels=["100","200","300","400","500","600","700","800","900","1000"])
  color_box(bp, 'red')
  plt.ylabel(r'$\bar{\alpha}_{P}$', fontsize=18)
  plt.xlabel('Number of samples', fontsize=18)
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.12,top=0.88)
  # SMC alpha_D
  plt.figure()
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=18)
  bp = plt.boxplot(MatrixPar4, showfliers=False, labels=["100","200","300","400","500","600","700","800","900","1000"])
  color_box(bp, 'red')
  plt.ylabel(r'$\bar{\alpha}_{D}$', fontsize=18)
  plt.xlabel('Number of samples', fontsize=18)
  plt.subplots_adjust(left=0.18,right=0.96,bottom=0.12,top=0.88)
  plt.show()


def PlotSamplesParameters():
  input1 = np.loadtxt("GPR_1000/samplesGPR.txt", dtype='f', delimiter=' ')
  input2 = np.loadtxt("../Output_AGPR/samples.txt", dtype='f', delimiter=' ')
  file1 = np.array(input1)
  file2 = np.array(input2)
  file1_Par1 = file1[:,0]
  file1_Par2 = file1[:,1]
  file2_Par1 = file2[:,0]
  file2_Par2 = file2[:,1]
  plt.scatter(file2_Par1,file2_Par2,color='blue',label='AGPR')
  plt.scatter(file1_Par1,file1_Par2,color='red',label='GPR')
  plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.14),ncol=6,fontsize=14)
  print(file2_Par1.shape)
  plt.xlabel(r'$\bar{\alpha}_{P}$', fontsize=18)
  plt.ylabel(r'$\bar{\alpha}_{D}$', fontsize=18)
  plt.show()

def Plot_MAP_response():
  input1 = np.loadtxt("MAPs/OutModel_100_MCMC_GPR.dat", dtype='f', delimiter='\t')
  input2 = np.loadtxt("MAPs/OutModel_100_SMC_GPR.dat", dtype='f', delimiter='\t')
  input3 = np.loadtxt("MAPs/OutModel_1000_MCMC_GPR.dat", dtype='f', delimiter='\t')
  input4 = np.loadtxt("MAPs/OutModel_1000_SMC_GPR.dat", dtype='f', delimiter='\t')
  QOI1 = np.array(input1)
  QOI2 = np.array(input2)
  QOI3 = np.array(input3)
  QOI4 = np.array(input4)
  Mean1 = np.mean(QOI1,axis=0)
  Mean2 = np.mean(QOI2,axis=0)
  Mean3 = np.mean(QOI3,axis=0)
  Mean4 = np.mean(QOI4,axis=0)  
  data = np.array([0.56829, 0.5681037, 0.5686626, 0.5691672, 0.6176569, 0.6802602, 0.7250002, 0.7600671, 0.7605562, 0.7590967, 0.7202258, 0.6275473, 0.538145, 0.4514444, 0.3496518, 0.2844012, 0.1959615, 0.1393592, 0.1014277, 0.06882953, 0.04892439, 0.0309368, 0.02281637, 0.01860866, 0.01670665, 0.01145865, 0.00999915, 0.01005349, 0.00989046, 0.00784871, 0.0077245, 0.00554301, 0.00337704, 0.01571294, 0.0155732, 0.0155111, 0.01561202, 0.0151928, 0.01500648, 0.01520833, 0.0155111, 0.01558873, 0.01829036, 0.07083246, 0.1857993, 0.2866292, 0.3938561, 0.4960524, 0.5698349, 0.6533526, 0.7113057, 0.7336408, 0.7633044, 0.7718518, 0.784506, 0.7873551, 0.7897307, 0.7894745, 0.7897229, 0.7902896, 0.7875647, 0.786703, 0.785818, 0.784638, 0.7825574, 0.7837685 ]) # [Live; Dead] -> parameter 0.2 and 0.15 with normal(0;0.03)
   
  time = np.linspace(0,96,33)
  plt.figure(figsize=(12,8))
  plt.errorbar(time,data[0:33],yerr=0.03,fmt='o',color='green',label='Live')
  plt.plot(time,Mean1[0:33],color='gray',label='MCMC (100)',ls='-.')
  plt.plot(time,Mean2[0:33],color='blue',label='SMC (100)',ls='-.')
  plt.plot(time,Mean3[0:33],color='gray',label='MCMC (1000)',ls=':')
  plt.plot(time,Mean4[0:33],color='blue',label='SMC (1000)',ls=':')
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
  plt.show()


# Read_Plot("GPR_100/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_100')
# Read_Plot("GPR_100/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_100')
# Read_Plot("GPR_200/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_200')
# Read_Plot("GPR_200/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_200')
# Read_Plot("GPR_300/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_300')
# Read_Plot("GPR_300/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_300')
# Read_Plot("GPR_400/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_400')
# Read_Plot("GPR_400/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_400')
# Read_Plot("GPR_500/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_500')
# Read_Plot("GPR_500/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_500')
# Read_Plot("GPR_600/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_600')
# Read_Plot("GPR_600/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_600')
# Read_Plot("GPR_700/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_700')
# Read_Plot("GPR_700/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_700')
# Read_Plot("GPR_800/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_800')
# Read_Plot("GPR_800/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_800')
# Read_Plot("GPR_900/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_900')
# Read_Plot("GPR_900/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_900')
# Read_Plot("GPR_1000/CalibSMC_GPR.dat",color="blue", name='Posterior_SMC_1000')
# Read_Plot("GPR_1000/CalibMCMC_GPR.dat",color="blue", name='Posterior_MCMC_1000')



# plt.show()

# Plot_boxplot()

# PlotSamplesParameters()

Plot_MAP_response() 

#Plot_STD("../Output_AGPR/diffenceSTD.txt")
