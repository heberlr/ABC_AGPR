
Run ABC methods will save parameters samples *.dat

Run AGPR will save Dictionary.pkl: {"GPR", "samples_parameters", "difference_std" }

Parallel version is for run just stochatic models, you need to define the number of task  equal  number of replicates + 1, example:

# Run the model hybrid with 10 replicates in parallel
``` mpiexec -n 11 py hybrid.py ```

version of libraries

Describe readme and look at Zenodo

Name: numpy
Version: 1.20.3

Name: scikit-learn
Version: 1.0.2

Name: lhsmdu
Version: 1.1

Name: pickle
Version: 4.0

Name: seaborn
Version: 0.10.0

Name: matplotlib
Version: 3.2.0

Name: scipy
Version: 1.7.1
