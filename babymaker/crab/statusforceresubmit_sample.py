#! /usr/bin/env python

# statusforceresubmit_sample.py: Checks status for forces re-submission of failed CRAB jobs

import os, sys

if len(sys.argv)<2: sys.exit('plese specify dataset')
dset        = sys.argv[1]

dsetnew = dset.split('/')[1]+'_'+dset.split('/')[2]
dsetnew = dsetnew.replace('SIM','')
dsetnew = dsetnew.replace('/AOD','')
dsetnew = dsetnew.replace('/USER','')
dsetnew = dsetnew.replace('/GEN-SIM-RAW','')
dsetnew = dsetnew.replace('/','_')


samples = ['el15_720elid']
for sam in samples:
    task = 'jobs/'+dsetnew+'_'+sam
    print('crab -status -c '+task+' >& tmp.log')
    os.system('crab -status -c '+task+' >& tmp.log')
    infile = open("tmp.log","r")
    lines = infile.readlines()

    failed = []
    check = 0
    for line in lines:
        myline = line.split()
        if len(myline)==7:
            if check:
                if myline[2]=="Done" or myline[1]=="Y":
                    #if myline[4]!="0" or myline[5]!="0": failed.append(myline[0]) 
                    if myline[4]!="0" or (myline[5]!="0" and myline[5]!="70500"): failed.append(myline[0]) 
            elif myline[0]=="-----": check=1
        elif len(myline)==0: check = 0

    if len(failed)>0:
        jlist = ""
        for job in failed:
            jlist = jlist+job+','
        jlist=jlist.rstrip(',')
        print('crab -forceResubmit '+jlist+' -c '+task)
        os.system('crab -forceResubmit '+jlist+' -c '+task)

    os.system('rm tmp.log') 
