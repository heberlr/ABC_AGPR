#! /usr/bin/env python3
#
import numpy as np
from mpi4py import MPI
import subprocess
from scipy.stats import multivariate_normal
import os

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
  
def ABC_SMC(Model, data, LowLimit, UpperLimit, FILE='CalibSMC.dat', tol = np.array([98,99,100]), NumAccept = 100, max_iterations=100000, var_trasition=0.2):
  #*****************************************************************************
  #
  ## Approximate  Bayesian Computation - Sequential Monte Carlo - Bayesian inference
  #  (Del Moral et al. - 2006 - Sequential Monte Carlo samplers)
  #
  #  Modified:
  #
  #  7 Set 2020
  #
  #  Author:
  #
  #    Heber L. Rocha
  #
  #  Input:
  #
  #    function Model(Par): model with output compatible to observational data.
  #
  #    real data[n]: contains the observational data
  #
  #    real LowLimit[numpar]: lower limit for parameters.
  #
  #    real UpperLimit[numpar]: upper limit for parameters.
  #
  #    string FILE: output file name.
  #
  #    array real tol: tolerance vector to evaluate the distance between observational and model data (must be in ascending order). The size of this vector is the number of populations.
  #    
  #    int NumAccept: the number of accepted parameters that it will generate a posterior distribution (for each population).
  #
  #    int max_iterations: the number max of execution of model for each population.
  #     
  #    real var_trasition: variance of the normal distribution for sampling. 
  
  Npar = UpperLimit.shape[0] # Number of parameters
  Nqoi = data.shape[0] # Number of quantity of interest
  if rank == 0:
      file = open(FILE,"w")    
      Npop = tol.shape[0]
      theta_star = np.zeros(Npar, dtype='d')
      theta_ant = np.zeros((NumAccept,Npar))
      weight_prev = np.ones(NumAccept)
      weight = np.ones(NumAccept)
      Kernel = np.zeros(NumAccept)
      #Loop populations
      for k in range(0, Npop):
        count = 0
        # Generate posterior of the population k
        for i in range(0, max_iterations):
            cond = True
            while(cond):
                if (k == 0):
                  for j in range(0, Npar):
                    theta_star[j] = np.random.uniform(LowLimit[j],UpperLimit[j])
                else:
                  index = np.random.choice(NumAccept,p=weight_prev)
                  theta_star = theta_ant[index,:] + np.random.normal(0, var_trasition*(UpperLimit-LowLimit))
                cond = [False for k in range(0,Npar) if theta_star[k]>UpperLimit[k] or theta_star[k]<LowLimit[k]]
            # Send and Receive data (MPI)
            QOI = np.zeros((size-1, Nqoi))
            for rankID in range(1,size):
                comm.Send(np.append(theta_star,0.0), dest=rankID, tag=rankID)
            for rankID in range(1,size):
                comm.Recv(QOI[rankID-1,:], source=rankID, tag=rankID+size)
            output_model = np.mean(QOI, axis=0)
            distance = np.sqrt(np.sum([(a - b)**2 for a, b in zip(output_model, data)]))
            print(str(count)+"/"+str(i)+" -- distance: "+str(distance)+" "+ str(theta_star)+"\n")
            if (distance < tol[k]):
                if (k == 0): weight[count] = 1.0
                else:
                  for j in range(0, NumAccept): Kernel[j] = np.linalg.norm(multivariate_normal.pdf(theta_ant[j,:], mean=theta_star, cov=var_trasition))
                  weight[count] = 1.0/np.sum(weight_prev*Kernel)
                # Add sample to population
                theta_ant[count,:] = theta_star
                count = count + 1
            if (count == NumAccept):
                break
        # Normalize weights 
        weight /= np.sum(weight)
        weight[np.isnan(weight)] = 0.0
        weight[-1] += np.abs(1.0-np.sum(weight))
        weight_prev = np.copy(weight)
        
      # Print last population
      for i in range(0, NumAccept):
        for j in range(0, Npar):
            file.write(str(theta_ant[i,j])+" ")
        file.write(str(weight[i])+"\n")  
      file.close()
      # Finished Threads
      for rankID in range(1,size):
        comm.Send(np.append(theta_star,1.0), dest=rankID, tag=rankID)
      return theta_ant
  else:
      Par = np.zeros(Npar+1, dtype='d')   
      while (Par[-1]==0.0):
        comm.Recv(Par, source=0,tag=rank)
        if (Par[-1] == 1.0): break
        OUT = Model(Par[:-1],Nqoi)
        comm.Send(OUT, dest=0,tag=rank+size)



def ABC_MCMC(Model, data, LowLimit, UpperLimit, FILE='CalibMCMC.dat', tol = 100, NumAccept = 100, max_iterations=100000, var_trasition=0.2):
  #*****************************************************************************
  #
  ## Markov chain Monte Carlo without likelihoods - Bayesian inference
  #  (Marjoram et al. - 2003 - Markov chain Monte Carlo without likelihoods)
  #
  #  Modified:
  #
  #  7 Set 2020
  #
  #  Author:
  #
  #    Heber L. Rocha
  #
  #  Input:
  #
  #    function Model(Par): model with output compatible to observational data.
  #
  #    real data[n]: contains the observational data
  #
  #    real LowLimit[numpar]: lower limit for parameters.
  #
  #    real UpperLimit[numpar]: upper limit for parameters.
  #
  #    string FILE: output file name.
  #
  #    real tol: tolerance between observational and model data.
  #    
  #    int NumAccept: the number of accepted parameters that it will generate a posterior distribution.
  #
  #    int max_iterations: the number max of execution of model.
  #     
  #    real var_trasition: variance of the normal distribution for sampling.
  
  Npar = UpperLimit.shape[0] # Number of parameters
  Nqoi = data.shape[0] # Number of quantity of interest
  if rank == 0:
      file = open(FILE,"w") 
      count = 0
      Npar = UpperLimit.shape[0]
      theta_star = np.zeros(Npar)
      theta = np.zeros((NumAccept,Npar))
      for j in range(0, Npar):
        theta_star[j] = np.random.uniform(LowLimit[j],UpperLimit[j])
      for i in range(0, max_iterations):
        # Send and Receive data (MPI)
        QOI = np.zeros((size-1, Nqoi))
        for rankID in range(1,size):
            comm.Send(np.append(theta_star,0.0), dest=rankID, tag=rankID)
        for rankID in range(1,size):
            comm.Recv(QOI[rankID-1,:], source=rankID, tag=rankID+size)
        output_model = np.mean(QOI, axis=0)  
        distance = np.sqrt(np.sum([(a - b)**2 for a, b in zip(output_model, data)]))
        print(str(count)+"/"+str(i)+" -- distance: "+str(distance)+" "+ str(theta_star)+"\n")
        if (distance < tol or count == 0):
            theta[count,:] = theta_star
            count = count + 1
            for j in range(0, Npar):
              file.write(str(theta_star[j])+" ")
            file.write(str(count)+" "+str(i)+" "+str(distance)+"\n")
        if (count == NumAccept):
            break
        cond = True
        while(cond):
          noise = np.random.normal(0, var_trasition*(UpperLimit-LowLimit))
          theta_star = theta[count-1,:] + noise
          cond = [False for k in range(0,Npar) if theta_star[k]>UpperLimit[k] or theta_star[k]<LowLimit[k]]
      file.close()
      # Finished Threads
      for rankID in range(1,size):
        comm.Send(np.append(theta_star,1.0), dest=rankID, tag=rankID)
      return theta
  else:
      Par = np.zeros(Npar+1, dtype='d')   
      while (Par[-1]==0.0):
        comm.Recv(Par, source=0,tag=rank)
        if (Par[-1] == 1.0): break
        OUT = Model(Par[:-1],Nqoi)
        comm.Send(OUT, dest=0,tag=rank+size)
  
def Calling_modelOpenMP(parameter,NumQOI):
    function_call = ['./build/main.exe', '{}'.format(rank)] #change of .exe to .out
    for ind in range(0,parameter.shape[0]):
        function_call.append('{}'.format(parameter[ind]))
    cache = subprocess.run(function_call,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(cache.stdout)
    # print(cache.stderr)
    # print(cache.returncode)
    #print(OUT)
    OUT = np.fromstring(cache.stdout, dtype='d', sep=' ')
    time = OUT[::3]
    Live = OUT[1::3]
    Dead = OUT[2::3]
    QOI = np.concatenate((Live, Dead), axis=None)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)
    if (QOI.shape[0] != NumQOI):
        print("Model output error! incompatibility of QoIs!")
        os._exit(0)
    return QOI     

def Calling_model(parameter,NumQOI):
    time = np.linspace(0,1,100)
    return parameter[0]*time + parameter[1]  

def  Test():
    data = 0.2*np.linspace(0,1,100) + 0.7
    UP = np.array([1.0,1.0])
    DW = np.array([0.0,0.0])
    ABC_MCMC(Calling_model,data,DW,UP,tol=1.0,NumAccept=100)    
    
data = np.array([0.56829, 0.5681037, 0.5686626, 0.5691672, 0.6176569, 0.6802602, 0.7250002, 0.7600671, 0.7605562, 0.7590967, 0.7202258, 0.6275473, 0.538145, 0.4514444, 0.3496518, 0.2844012, 0.1959615, 0.1393592, 0.1014277, 0.06882953, 0.04892439, 0.0309368, 0.02281637, 0.01860866, 0.01670665, 0.01145865, 0.00999915, 0.01005349, 0.00989046, 0.00784871, 0.0077245, 0.00554301, 0.00337704, 0.01571294, 0.0155732, 0.0155111, 0.01561202, 0.0151928, 0.01500648, 0.01520833, 0.0155111, 0.01558873, 0.01829036, 0.07083246, 0.1857993, 0.2866292, 0.3938561, 0.4960524, 0.5698349, 0.6533526, 0.7113057, 0.7336408, 0.7633044, 0.7718518, 0.784506, 0.7873551, 0.7897307, 0.7894745, 0.7897229, 0.7902896, 0.7875647, 0.786703, 0.785818, 0.784638, 0.7825574, 0.7837685 ]) # [Live; Dead] -> parameter 0.2 and 0.15 with normal(0;0.03)
UP = np.array([0.5,0.5])
DW = np.array([0.0,0.0])
ABC_SMC(Calling_modelOpenMP,data,DW,UP,tol=np.array([1.0,0.5,0.2]),NumAccept=1000)
ABC_MCMC(Calling_modelOpenMP,data,DW,UP,tol=0.2,NumAccept=1000,var_trasition=0.1)


