#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:16:50 2014

@author: patricio
"""

import sys

deps= {"bitarray":"bitarray", "Bio":"biopython", "biom":"biom-format"}

ok=True

for dep in deps.iterkeys():
    try:
        __import__(dep)
    except ImportError:
        print("Required module {} not installed. Please install.".format(deps[dep]))
        ok=False
    else:
        print("{} imported succesfully".format(deps[dep]))

if ok:
    print("All required modules imported successfully.")
else:
    print("Some python modules were not found. Please install them before running the pipeline")
    sys.exit(1)

sys.exit(0)