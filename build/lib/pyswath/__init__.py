######!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import  division
# This to be sure that the result of the division of integers is a real, not an integer

__version__ = "0.1.0"

# Import modules
import sys
import os
#import copy
import numpy as np

# Import all the functions
__all__ = ['raster_tools', 'profiles', 'swath', 'plotgraph', 'checks']

from raster_tools import *
from profiles import *
from swath import *
from plotgraph import *
from checks import *
