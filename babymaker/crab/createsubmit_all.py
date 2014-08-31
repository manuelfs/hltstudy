#! /usr/bin/env python

# creatsubmit_all.py: Creates and submits CRAB jobs based on a file with datasets

import os, sys

if len(sys.argv)<2: sys.exit('Please specify the filename with the datasets')

infile = open(sys.argv[1],"r")
lines = infile.readlines() 

for line in lines:
    if line.find('#')==-1 :
        mylineWithXs = line.split()
        if len(mylineWithXs)>0:
            myline = mylineWithXs[0]
            if myline!='':
                print('-------------------------------------------------------------------------------------------------------------------')
                print('./CrabAuto/createsubmit_sample.py '+line)
                os.system('./CrabAuto/createsubmit_sample.py '+line)
                print('-------------------------------------------------------------------------------------------------------------------\n')
