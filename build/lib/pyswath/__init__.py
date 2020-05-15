######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later


from __future__ import  division
# This to be sure that the result of the division of 2 integers is a real, not an integer
from __future__ import absolute_import
from __future__ import print_function

__version__ = "2.0.4"

# Import modules
import sys
import os
#import copy
import numpy as np

# Import all the functions
__all__ = ['raster_tools', 'profiles', 'swath', 'plotgraph', 'checks']

from .raster_tools import *
from .profiles import *
from .swath import *
from .plotgraph import *
from .checks import *
