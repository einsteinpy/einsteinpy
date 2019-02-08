import numpy as np
import pytest
from astropy import units as u
from numpy.testing import assert_allclose

from einsteinpy import constant
from einsteinpy.metric import Schwarzschild

_c = constant.c.value

