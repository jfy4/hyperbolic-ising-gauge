#!/home/jfunmuth/miniconda3/bin/python

import os
import sys
import time
import numpy as np


L = 18
#betas = np.round(np.asarray([0.1 + 0.1*i for i in range(28)] + [3.2 + 0.4*i for i in range(8)]), 1)
# betas = np.round([0.4 + 0.4*i for i in range(15)], 1)
betas = np.round([2.0, 3.6, 5.6], 1)
# betas = np.round(np.asarray([0.1, 0.2, 0.3, 0.5, 0.6, 0.7]), 1)
for b in betas:
    os.system("mkdir hb_b" + str(b).replace('.', 'p'))
    os.system("mkdir hb_b" + str(b).replace('.', 'p') + "/L" + str(L))	
    with open("job.submit", 'r') as file:
        text = file.readlines()
        
        text[1] = "arguments               = " + str(b) + "\n"
        
    with open("job.submit", 'w') as file:
        file.writelines(text)
    os.system("condor_submit job.submit -batch-name " + "hbrunL" + str(L) + "b" + str(b))
    time.sleep(2.)    

