#!/home/jfunmuth/miniconda3/bin/python

import os

for i in range(3971082-3971068):
    os.system("condor_rm " + str(3971068 + i) + ".0")
 
