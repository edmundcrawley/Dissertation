"""
 Runs code for all three chapters of the dissertation
"""

# Requrements:
#     numdifftools
#     numpy
#     matplotlib
#     scipy

import sys 
import os

Rscript_path = "C:/Program Files/R/R-3.4.3/bin/x64/Rscript"
DissertationPath = "C:/Users/edmun/OneDrive/Documents/Research/Dissertation/"
sys.path.append(DissertationPath + 'Chapter1/')

#Chapter 1
os.chdir(DissertationPath + "Chapter1/")
import do_Chapter1
