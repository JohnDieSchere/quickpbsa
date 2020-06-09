#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
submodule utilities
===================

Provides some general utilities

Author: Johan Hummert

If you find this software useful, please reference: REF_TO_PAPER
"""

# relative imports 
# (try ... except to skip parts if tifffile / matplotlib is not installed)
try:
    from .plot_selection import *
except:
    pass
try:
    from .stacks import *
except:
    pass



