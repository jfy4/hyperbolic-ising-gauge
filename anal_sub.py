#!/home/jfunmuth/miniconda3/bin/python

import os
import sys
import time
import numpy as np


# betas = np.round(np.asarray([0.1 + 0.1*i for i in range(28)] + [3.2 + 0.4*i for i in range(8)]), 1)
betas = np.round(np.asarray([0.4 + 0.4*i for i in range(15)]), 1)
# betas = np.round([2.0, 3.6, 5.6], 1)
# masses = [0.004]
for b in betas[1:]:
    with open("job.submit.anal", 'r') as file:
        text = file.readlines()
        
        text[1] = "arguments               = " + str(b) + "\n"
        
    with open("job.submit.anal", 'w') as file:
        file.writelines(text)
    os.system("condor_submit job.submit.anal -batch-name " + "analL18" + str(b).replace('.', 'p'))
    time.sleep(2.)    
