#! /usr/bin/env python

# statusforceresubmit_all.py: Checks status for forces re-submission of failed CRAB jobs from file

import os, sys

if len(sys.argv)<2: sys.exit('plese specify the filename')

infile = open(sys.argv[1],"r")
lines = infile.readlines() 

for line in lines:
    if line.find('#')==-1 :
        mylineWithXs = line.split()
        if len(mylineWithXs)>0:
            myline = mylineWithXs[0]
            if myline!='':
                print('-------------------------------------------------------------------------------------------------------------------')
                print('./CrabAuto/statusforceresubmit_sample.py '+line)
                os.system('./CrabAuto/statusforceresubmit_sample.py '+line)
                print('-------------------------------------------------------------------------------------------------------------------\n')
