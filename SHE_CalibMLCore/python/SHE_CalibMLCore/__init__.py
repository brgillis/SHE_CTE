"""
This code is based on Tenbilac by Malte Tewes and Thibault Kuntzer
https://github.com/mtewes/tenbilac

It provides a numpy implementation of a feedforward neural network designed to solve an inverse regression problem 
(aka "calibration problem" of regression). It performs a supervised machine learning regression from a noisy feature 
space to a well known explanatory space, and can be trained to minimize bias instead of error in this explanatory space.
"""

from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)


# We add "aliases" to sys.modules using the "tenbilac" names, so that SHE_CalibMLCore can unpickle tenbilac objects
# For this we first import the required modules
from . import act
from . import com
from . import data
from . import err
from . import layer
from . import net
from . import train

import sys

sys.modules["tenbilac.train"] = sys.modules["SHE_CalibMLCore.train"]
sys.modules["tenbilac.layer"] = sys.modules["SHE_CalibMLCore.layer"]
sys.modules["tenbilac.err"] = sys.modules["SHE_CalibMLCore.err"]
sys.modules["tenbilac.net"] = sys.modules["SHE_CalibMLCore.net"]
sys.modules["tenbilac.act"] = sys.modules["SHE_CalibMLCore.act"]
sys.modules["tenbilac.data"] = sys.modules["SHE_CalibMLCore.data"]
