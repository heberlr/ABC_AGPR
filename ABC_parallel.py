def ABC_SMC(Model, data, LowLimit, UpperLimit, FILE='CalibSMC.dat', tol = ([98,99,100]), NumAccept = 1000, max_iterations=100000, var_trasition=0.2):
  #*****************************************************************************
  #
  ## Approximate  Bayesian Computation - Sequential Monte Carlo - Bayesian inference
  #  (Del Moral et al. - 2006 - Sequential Monte Carlo samplers)
  #
  #  Modified:
  #
  #  29 May 2020
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
  
  from mpi4py import MPI
  import numpy as np
  import subprocess
  from scipy.stats import multivariate_normal
  
  comm = MPI.COMM_WORLD
  size = comm.Get_size()
  rank = comm.Get_rank()
  
  file = open(FILE,"w") 
  Npar = UpperLimit.shape[0]
  Npop = tol.shape[0]
  theta_star = np.zeros(Npar)
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
        output_model = Model(theta_star)
        distance = np.sqrt(np.sum([(a - b)**2 for a, b in zip(output_model, data)]))
        print(str(count)+"/"+str(i)+" -- distance: "+str(distance)+" "+ str(theta_star)+"\n")
        if (distance < tol[k]):
            if (k == 0): weight[count] = 1.0
            else:
              for j in range(0, NumAccept): Kernel[j] = np.linalg.norm(multivariate_normal.pdf(theta_ant[j,:], mean=theta_star, cov=var_trasition))
              weight[count] = 1.0/np.sum(weight_prev*Kernel)
              weight[np.isnan(weight)] = 0.0
            # Add sample to population
            theta_ant[count,:] = theta_star
            count = count + 1
        if (count == NumAccept):
            break
    # Normalize weights 
    weight /= np.sum(weight)
    weight[-1] += np.abs(1.0-np.sum(weight))
    weight_prev = np.copy(weight)

  # Print last population
  for i in range(0, NumAccept):
    for j in range(0, Npar):
      file.write(str(theta_ant[i,j])+" ")
    file.write(str(weight[i])+"\n")  
  file.close()
  return theta_ant 

def SampleTheta(NumProcess, LowLimit, UpperLimit, NumAccept, weight_prev, theta_ant, var_trasition):
  Npar = UpperLimit.shape[0]
  theta = np.zeros((NumProcess,Npar))  
  for i in range(0, NumProcess):
    cond = True
    while(cond):
      if (k == 0):
        for j in range(0, Npar):
          theta_star[j] = np.random.uniform(LowLimit[j],UpperLimit[j])
      else:
        index = np.random.choice(NumAccept,p=weight_prev)
        theta_star = theta_ant[index,:] + np.random.normal(0, var_trasition*(UpperLimit-LowLimit))
      cond = [False for k in range(0,Npar) if theta_star[k]>UpperLimit[k] or theta_star[k]<LowLimit[k]]
    theta[i,:] = theta_star
  return theta
