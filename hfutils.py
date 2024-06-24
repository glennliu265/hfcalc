#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Heat Flux Feedback Calculation Utilities
Functions are largely copied from amv.proc. Copy date is indicated in documentation
for each function.

Functions
    
    addstrtoext : Append string to end of file, before the extension

Created on Fri Jun 21 15:47:29 2024

@author: gliu

"""

def addstrtoext(name,addstr,adjust=0):
    """
    Copied from amv.proc on 2024.06.24
    
    Add [addstr] to the end of a string with an extension [name.ext]
    Result should be "name+addstr+.ext"
    -4: 3 letter extension. -3: 2 letter extension
    """
    return name[:-(4+adjust)] + addstr + name[-(4+adjust):]