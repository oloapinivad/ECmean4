"""ECmean4 module"""

__version__ = '0.1.11'

# functions to be accessible everywhere
from ecmean.libs.diagnostic import Diagnostic
from ecmean.libs.support import Supporter
from ecmean.libs.units import UnitsHandler
from ecmean.global_mean import global_mean
from ecmean.performance_indices import performance_indices

__all__ = ["global_mean", "performance_indices", "Diagnostic", "Supporter", "UnitsHandler"]
