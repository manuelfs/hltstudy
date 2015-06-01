#! /usr/bin/env python

# creatsubmit_sample.py: Creates and submits individual CRAB jobs

import os, sys

if len(sys.argv)<2: sys.exit('plese specify the dataset')

dset = sys.argv[1]

#if len(sys.argv)==3:
#    print 'using gtag='+sys.argv[2]
#    gtag = sys.argv[2]

dsetnew = dset.split('/')[1]+'_'+dset.split('/')[2]
dsetnew = dsetnew.replace('SIM','')
dsetnew = dsetnew.replace('/AOD','')
dsetnew = dsetnew.replace('/USER','')
dsetnew = dsetnew.replace('/GEN-SIM-RAW','')
dsetnew = dsetnew.replace('/','_')


samples = ['mu15_743', 'el15_743']
for sam in samples:
    task = 'jobs/'+dsetnew+'_'+sam
    os.system('rm -rf '+task)
    os.system('./makecfg3.py -d '+dset+' -s '+sam)  
    #print('crab -create -cfg cfg/'+dsetnew+'_'+sam+'.cfg\n')
    #print('crab -submit -c '+task) 
    #os.system('crab -create -cfg cfg/'+dsetnew+'_'+sam+'.cfg')
    #os.system('crab submit -c '+task) 
    print('crab submit -c cfg/'+dsetnew+'_'+sam+'.py') 
    os.system('crab submit -c cfg/'+dsetnew+'_'+sam+'.py') 

#if isdata: os.system('./makecfg.py -CMS2cfg ../src/CMS2/NtupleMaker/test/Slim_Data53X_SDFilter_crab_cfg.py -d  '+dset+' -t '+tag+' -evtsPerJob '+nevts+' -gtag '+gtag)
#else: os.system('./makecfg.py -CMS2cfg ../src/CMS2/NtupleMaker/test/Slim_MCProduction2012_SingleOrDiLepFilter_cfg.py -d  '+dset+' -t '+tag+' -evtsPerJob '+nevts+' -gtag '+gtag)
#if isdata: os.system('perl -p -i -e "s/events/lumis/g" '+dsetnew+'.cfg')
#os.system('crab -create -cfg '+dsetnew+'.cfg')
#os.system('crab -submit -c '+dsetnew) 
