#Just a quick script to check if any of the stored generic vibrations has
#an imaginary mode.

import pickle
import os
import numpy as np

os.chdir('generic-vibrations')
A = os.listdir('.')
for file in A:
    B = pickle.load(open(file))
    vibs = B['vibrations']
    count =0
    for vib in vibs:
        if np.imag(vib)>0:
            count=count+1
    if count>0:
        print file, ' has ', str(count), ' imaginary modes'
    
