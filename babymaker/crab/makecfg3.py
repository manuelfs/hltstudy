#! /usr/bin/env python

# makecfg.py: Creates configuration file for CRAB job

import string
import commands, os, re
import sys 
                      

dataSet = ''
outNtupleName = 'ntuple_hlt'
sample   = '';

def makeCrabConfig():
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2]
    jobname = dataSet[1:].replace('/','__') + '__' + sample
    jobname = jobname.replace(':','___')
    fullFileName = 'cfg/'+outFileName + '_' + sample +'.py'
    outFile = open(fullFileName, 'w')
    print 'Writing CRAB config file: ' + fullFileName


    outFile.write("from WMCore.Configuration import Configuration                  \n")
    outFile.write("config = Configuration()                                        \n")
    outFile.write("                                                                \n")
    outFile.write("config.section_('General')                                      \n")
    outFile.write("config.General.requestName = '"+jobname+"'                      \n")
    outFile.write("config.General.workArea = 'out_crab'                            \n")
    outFile.write("                                                                \n")
    outFile.write("config.section_('JobType')                                      \n")
    outFile.write("config.JobType.pluginName = 'Analysis'                          \n")
    outFile.write("config.JobType.psetName = '../python/cfg_" + sample + ".py'     \n")
    outFile.write("                                                                \n")
    outFile.write("config.section_('Data')                                         \n")
    outFile.write("config.Data.inputDataset = '"+dataSet+"'                        \n")
    outFile.write("config.Data.inputDBS = 'global'                                 \n")
    outFile.write("config.Data.splitting = 'FileBased'                             \n")
    outFile.write("config.Data.unitsPerJob = 10                                    \n")
    outFile.write("config.Data.publication = True                                  \n")
    outFile.write("config.Data.publishDBS = 'phys03'                               \n")
    outFile.write("                                                                \n")
    outFile.write("config.section_('Site')                                         \n")
    outFile.write("config.Site.storageSite = 'T2_US_UCSD'                          \n")


for i in range(0, len(sys.argv)):
    if sys.argv[i] == '-d':
        dataSet = sys.argv[i+1]
    if sys.argv[i] == '-s':
        sample = sys.argv[i+1]

makeCrabConfig()

