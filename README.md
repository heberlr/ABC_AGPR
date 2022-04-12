The folder ```./src``` contain all methods code with version serial and parallels. The version in parallel is for run just stochastic models, where you need to define number of replicates for calibrate the mean of simulations.

The folders ```./Hybrid``` and ```./Verhulst``` contain scripts (```hybrid.py``` and ```verhulst.py```) to generate the AGPR and executing the methods ABC MCMC and ABC SMC with and without AGPR. The ABC methods will save sample parameters that generate the posterior distribution in ```.dat``` files inside the ```./Calibration``` folder. The AGPR will save a file ```Dictionary.pkl``` in ```./AGPR``` composed by ```{"GPR","samples_parameters","difference_std"}```.
- "GPR" - model trained to generate predictions
- "samples_parameters" - samples of parameters used for training
- "difference_std" - difference between max(std) and mean(std) of elements in parametric hypercube.

In case of the hybrid model, the script ```hybrid.py``` needs to be executed in parallel using more than one task. The number of tasks must be equal to the number of replicates + 1 of the hybrid model. For example, if you want to define 10 replicates to run the hybrid model, you can define 11 tasks (``` mpiexec -n 11 py hybrid.py ```).

To visualize the results you can use the ```plot.py``` scripts.


Required packages:  (versions used)

* NumPy (version: 1.20.3)
* scikit-learn (version: 1.0.2)
* lhsmdu (version: 1.1)
* pickle (version: 4.0)
* seaborn (version: 0.10.0)
* matplotlib (version: 3.2.0)
* scipy (version: 1.7.1)
* mpi4py (version: 3.0.3)
