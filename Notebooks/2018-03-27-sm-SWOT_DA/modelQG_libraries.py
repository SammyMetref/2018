
# Model specific libraries
from importlib.machinery import SourceFileLoader 


modgrid = SourceFileLoader("modgrid", "/home/metrefs/Boost-SWOT/2018/Notebooks/Models/modgrid.py").load_module() 

moddyn = SourceFileLoader("moddyn", "/home/metrefs/Boost-SWOT/2018/Notebooks/Models/moddyn.py").load_module() 

modelliptic = SourceFileLoader("modelliptic", "/home/metrefs/Boost-SWOT/2018/Notebooks/Models/modelliptic.py").load_module() 

qgsw = SourceFileLoader("qgsw", "/home/metrefs/Boost-SWOT/2018/Notebooks/Models/qgsw.py").load_module()  



