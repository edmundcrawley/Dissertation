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
from time import clock
mystr = lambda number : "{:.4f}".format(number)

global Rscript_path 
Rscript_path = "C:/Program Files/R/R-3.4.3/bin/x64/Rscript"

DissertationPath = "C:/Users/edmun/OneDrive/Documents/Research/Dissertation/"
sys.path.append(DissertationPath + 'Chapter1/')

#Chapter 1
start_time = clock()
os.chdir(DissertationPath + "Chapter1/")
import do_Chapter1
end_time = clock()
print('Code for Chapter 1 took ' + mystr((end_time-start_time)/60.0) + ' minutes.')


sys.path.append(DissertationPath + 'Chapter2/Code/')

#Chapter 1
start_time = clock()
os.chdir(DissertationPath + "Chapter2/Code/")
import do_Chapter2
end_time = clock()
print('Code for Chapter 2 took ' + mystr((end_time-start_time)/60.0) + ' minutes.')

