import sys
sys.path.append('../src')
from plot import PlotPosterior


PlotPosterior("CalibMCMC.dat",2,color="gray")