import numpy as np
import sys

# Fixed seed - parameters: growth rate = 0.05 and carrying capacity = 10000.0
rng = np.random.RandomState(1234)
Par = np.array([0.05,10000.0])
t = np.linspace(24,288,12) #1-12 days
OutputModel = np.zeros(t.size)
InitialCond = 1000.0
std = InitialCond*0.2
OutputModel = (InitialCond*Par[1]*np.exp(Par[0]*t))/(InitialCond*(np.exp(Par[0]*t)-1) + Par[1]) + rng.normal(0, std, t.size)

np.savetxt("modeldata.dat", np.column_stack( (t, OutputModel, np.ones(t.size) * std ) ), delimiter = '\t')
np.savetxt("trueparams.dat", Par)
np.savetxt("ic.dat", [InitialCond])

