import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, DotProduct, Matern, ExpSineSquared


# Obtain model 
data = np.loadtxt("modeldata.dat", delimiter = '\t')
t, qoi, qoi_std = data[:,0], data[:, 1], data[:, 2]
theta_st = np.loadtxt("trueparams.dat").reshape(1,-1)
InitialCond = np.loadtxt("ic.dat")
UpperLimit = [0.11, 15000.0]
LowLimit = [0.0, 5000.0]
samples = 200

# Range of parameter, used in surrogate model construction
theta = np.array([np.random.uniform(LowLimit[0],UpperLimit[0],samples),np.random.uniform(LowLimit[1],UpperLimit[1],samples)]).T
# Model output associated with sampled parameters
f = np.zeros( (theta.shape[0],t.size) )
for i in range(0, theta.shape[0]):
  f[i,:] = (InitialCond*theta[i,1]*np.exp(theta[i,0]*t))/(InitialCond*(np.exp(theta[i,0]*t)-1) + theta[i,1])



# Creating GP kernels
hyp = np.array([0.1*(UpperLimit[0]-LowLimit[0]),0.1*(UpperLimit[1]-LowLimit[1])])
kernelRBF = ConstantKernel( np.random.rand()* 2, (1e-5, 1e+5) ) * RBF(hyp)
kernelMatern = ConstantKernel( np.random.rand()* 2, (1e-5, 1e+7) ) * Matern(hyp)
kernelDot = DotProduct( sigma_0 = 2.4, sigma_0_bounds=(1e-5, 1e9) )  
kernelExpSin = ConstantKernel( np.random.rand()* 2, (1e-5, 1e+7) ) * ExpSineSquared()


# Creating GP Regressor (RBF)
gpr = GaussianProcessRegressor(kernel=kernelRBF, n_restarts_optimizer=33, alpha=2 * 10**-3) 
gpr.fit(theta,f)
output = gpr.predict(theta_st, return_std=True)
print(output)
# Saving QoI, QoI Std and hyperparameters
np.savetxt("sqrexp_qoi.txt", np.squeeze(output[0]))
np.savetxt("sqrexp_qoi_std.txt", output[1])
with open('sqrexp_hypers.txt', "w") as file:
    file.write(str(gpr.kernel))


# Creating GP Regressor (Matern)
gpr = GaussianProcessRegressor(kernel=kernelMatern, n_restarts_optimizer=33, alpha=2 * 10**-3) 
gpr.fit(theta,f)
output = gpr.predict(theta_st, return_std=True)
# Saving QoI, QoI Std and hyperparameters
np.savetxt("Matern_qoi.txt", np.squeeze(output[0]))
np.savetxt("Matern_qoi_std.txt", output[1])
with open('Matern_hypers.txt', "w") as file:
    file.write(str(gpr.kernel))

# Creating GP Regressor (Dot)
gpr = GaussianProcessRegressor(kernel=kernelDot, n_restarts_optimizer=33, alpha=2 * 10**-3) 
gpr.fit(theta,f)
output = gpr.predict(theta_st, return_std=True)
# Saving QoI, QoI Std and hyperparameters
np.savetxt("Dot_qoi.txt", np.squeeze(output[0]))
np.savetxt("Dot_qoi_std.txt", output[1])
with open('Dot_hypers.txt', "w") as file:
    file.write(str(gpr.kernel))

# Creating GP Regressor (ExpSinSqred)
gpr = GaussianProcessRegressor(kernel=kernelExpSin, n_restarts_optimizer=33, alpha=2 * 10**-3) 
gpr.fit(theta,f)
output = gpr.predict(theta_st, return_std=True)
# Saving QoI, QoI Std and hyperparameters
np.savetxt("Expsinsquared_qoi.txt", np.squeeze(output[0]))
np.savetxt("Expsinsquared_qoi_std.txt", output[1])
with open('Expsinsquared_hypers.txt', "w") as file:
    file.write(str(gpr.kernel))