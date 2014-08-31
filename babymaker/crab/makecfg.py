#! /usr/bin/env python

# makecfg.py: Creates configuration file for CRAB job

import string
import commands, os, re
import sys 
                      

cmsswSkelFile = ''
dataSet = ''
numEvtsTotal = -1
numEvtsPerJob = 10000
outNtupleName = 'ntuple_hlt'
storageElement = 'T2_US_UCSD'
tag = 'V01-02-06'
mode = 'remoteGlidein'
dbs_url = ''
dbs_url_pub = 'http://ming.ucsd.edu:8080/DBS2/servlet/DBSServlet';
report_every = 1000;
global_tag_flag = '';
sParms = [];

fastSim = False;
MCatNLO = False;
isData  = False;
sample   = '';

def makeCrabConfig():
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2]
    fullFileName = 'cfg/'+outFileName + '_' + sample +'.cfg'
    outFile = open(fullFileName, 'w')
    print 'Writing CRAB config file: ' + fullFileName
    outFile.write('[CRAB]\n')
    outFile.write('jobtype                  = cmssw\n')
    outFile.write('scheduler                = ' + mode + '\n')
    outFile.write('use_server               = ' + '0' + '\n')
    outFile.write('\n[CMSSW]\n')
    outFile.write('datasetpath              = ' + dataSet + '\n')
    outFile.write('pset                     = ../python/run_' + sample + '.py \n')
    outFile.write('total_number_of_events   = ' + str(numEvtsTotal) + '\n')
    outFile.write('events_per_job           = ' + str(numEvtsPerJob) + '\n')
    outFile.write('output_file              = ' + outNtupleName + '_' + sample +'.root\n')
#    outFile.write('use_dbs3                 = 1\n')
    if dbs_url != '' :
        outFile.write('dbs_url                 = ' + dbs_url + '\n')
    if isData == True :
        outFile.write('total_number_of_lumis   = -1\n')
        outFile.write('lumis_per_job           = 1\n')
    outFile.write('\n')
    outFile.write('[USER]\n')
    outFile.write('return_data             = 0\n')
    outFile.write('copy_data               = 1\n')
    outFile.write('storage_element         = ' + storageElement + '\n')
    outFile.write('ui_working_dir          = jobs/' + outFileName + '_'+sample+'\n')
    outFile.write('user_remote_dir         = HLT/' + outFileName + '_'+sample+'\n')
    outFile.write('publish_data            = 0\n')
#    outFile.write('publish_data_name       = CMS2_' + tag + '\n')
#    outFile.write('dbs_url_for_publication = ' + dbs_url_pub + '\n\n')
    outFile.write('\n')
    outFile.write('[GRID]\n')
    # outFile.write('maxtarballsize = 20\n')
    outFile.write('data_location_override = T2_US\n')
    
	 
#    outFile.write('##here are some default sites that we \n')
#    outFile.write('##run at. Comment/Uncomment at will\n')
#    outFile.write('##UCSD \n')
#    outFile.write('#SE_white_list = T2_US_UCSD\n')
#    outFile.write('##WISC\n')
#    outFile.write('#SE_white_list = T2_US_Wisconsin\n')
#    outFile.write('##DESY\n')
#    outFile.write('#SE_white_list = T2_DE_DESY\n')
#    outFile.write('##Purdue\n')
#    outFile.write('#SE_white_list = T2_US_Purdue\n')
#    outFile.write('##MIT\n')
#    outFile.write('#SE_white_list = T2_US_MIT\n')
#    outFile.write('##Nebraska\n')
#    outFile.write('#SE_white_list = T2_US_Nebraska\n')
#    outFile.write('##IFCA\n')
#    outFile.write('#SE_white_list = T2_ES_IFCA\n')
#    outFile.write('##Lyon\n')
#    outFile.write('#SE_white_list = T2_FR_CCIN2P3\n')
#    outFile.write('##CIEMAT\n')
#    outFile.write('#SE_white_list = T2_ES_CIEMAT\n')
#    outFile.write('##IIHE\n')
#    outFile.write('#SE_white_list = T2_BE_IIHE\n')
#    outFile.write('##Aachen\n')
#    outFile.write('#SE_white_list = T2_DE_RWTH\n')
#
#def makeCMSSWConfig(cmsswSkelFile):
#    foundOutNtupleFile = False
#    foundreportEvery   = False
#    foundcmsPath       = False
#    inFile = open(cmsswSkelFile, 'r').read().split('\n')
#    nlines = 0
#    iline  = 0
#    for i in inFile:
#        nlines += 1
#        if i.find(outNtupleName) != -1:
#            foundOutNtupleFile = True
#        if i.find('reportEvery') != -1:
#            foundOutNtupleFile = True
#    if foundOutNtupleFile == False:
#        print 'The root file you are outputting is not named ntuple.root as it should be for a CMS2 job.'
#        print 'Please check the name of the output root file in your PoolOutputModule, and try again'
#        print 'Exiting!'
#        sys.exit()
#    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2] + '_cfg.py'
#    print 'Writing CMS2 CMSSW python config file : ' + outFileName
#    outFile = open(outFileName, 'w')
#    outFile.write( 'import sys, os' + '\n' + 'sys.path.append( os.getenv("CMSSW_BASE") + "/src/CMS2/NtupleMaker/test" )' + '\n' )
#    for i in inFile:
#        iline += 1
#        if i.find('reportEvery') != -1:
#            outFile.write('process.MessageLogger.cerr.FwkReport.reportEvery = ' + str(report_every) + '\n'); continue
#
#        if i.find('globaltag') != -1:
#            outFile.write('process.GlobalTag.globaltag = "' + global_tag + '"\n'); continue
#
#        if (i.find('cms.Path') != -1 and foundcmsPath == False ):
#            foundcmsPath = True            
#            outFile.write('process.eventMaker.datasetName                   = cms.string(\"' + dataSet+'\")\n')
#            outFile.write('process.eventMaker.CMS2tag                       = cms.string(\"' + tag+'\")\n')
#
#        outFile.write(i)
#        if iline < nlines:
#            outFile.write('\n')
#
#    # if foundcmsPath == True:
#    #   outFile.write('process.eventMaker.datasetName                   = cms.string(\"' + dataSet+'\")\n')
#    #   outFile.write('process.eventMaker.CMS2tag                       = cms.string(\"' + tag+'\")\n')
#      
#    if len(sParms) > 0:
#        outFile.write('process.sParmMaker.vsparms = cms.untracked.vstring(\n')
#        for sParm in sParms:
#            if sParm != sParms[-1]:  #assumes the list is populated with unique entries
#                sParm = '\"%s\",'%sParm
#            else:
#                sParm = '\"%s\"'%sParm
#            outFile.write('%s\n'%sParm)
#        outFile.write(') # list of sparm parameters, be sure it is the same size as the number of parameter in the files\n')
#        outFile.write('process.cms2WithEverything.replace( process.eventmakers, process.eventmakerswsparm ) #adds the sparm producer in to the sequence\n')
#
#    if fastSim:
#        outFile.write('process.cms2WithEverything.remove( process.cms2HBHENoiseFilterResultProducer ) #need to remove this produce for fastsim\n')
#
#    if MCatNLO:
#        outFile.write('process.cms2WithEverything.remove(process.CMS2FlavorHistorySequence)\n')
#
#
#    outFile.close()
#
#
#
#       
#
#
#if len(sys.argv) < 5 :
#    print 'Usage: makecfg.py [OPTIONS]'
#    print '\nWhere the required options are: '
#    print '\t-CMS2cfg\tname of the skeleton CMS2 config file '
#    print '\t-d\t\tname of dataset'
#    print '\t-t\t\tCMS2 tag, will be added to publish_data_name'
#    print '\nOptional arguments:'
#    print '\t-isData\t\tFlag to specify if you are running on data.'
#    print '\t-strElem\tpreferred storage element. Default is T2_US_UCSD if left unspecified'
#    print '\t-nEvts\t\tNumber of events you want to run on. Default is -1'
#    print '\t-evtsPerJob\tNumber of events per job. Default is 20000'
#    #print '\t-n\t\tName of output Ntuple file. Default is ntuple.root'
#    print '\t-m\t\tsubmission mode (possible: condor_g, condor, glite). Default is glidein'
#    print '\t-dbs\t\tdbs url'
#    print '\t-re\t\tMessage Logger modulus for error reporting. Default is 1000'
#    print '\t-gtag\t\tglobal tag. Default is MC_31X_V3::All'
#    print '\t-sParms\t\tComma seperated, ordered list of Susy Parameter names.'
#    print '\t-fastSim\t\tUse a subset of the sequence that is compatible with FastSim. Default is to not use it.'
#    print '\t-MCatNLO\t\tUse a subset of the sequence that is compatible with MC@NLO samples. Default is to not use it.'
#    sys.exit()
#
#
for i in range(0, len(sys.argv)):
    if sys.argv[i] == '-CMS2cfg':
        cmsswSkelFile = sys.argv[i+1]
    if sys.argv[i] == '-d':
        dataSet = sys.argv[i+1]
    if sys.argv[i] == '-nEvts':
        numEvtsTotal = sys.argv[i+1]
    if sys.argv[i] == '-evtsPerJob':
        numEvtsPerJob = sys.argv[i+1]
    if sys.argv[i] == '-strElem':
        storageElement = sys.argv[i+1]
    #if sys.argv[i] == '-n':
    #    outNtupleName  = sys.argv[i+1]
    if sys.argv[i] == '-t':
        tag  = str(sys.argv[i+1])
    if sys.argv[i] == '-m':
        mode  = str(sys.argv[i+1])
    if sys.argv[i] == '-dbs':
        dbs_url = str(sys.argv[i+1])
    if sys.argv[i] == '-re':
        report_every = str(sys.argv[i+1])
    if sys.argv[i] == '-gtag':
        global_tag_flag = str(sys.argv[i+1])
    if sys.argv[i] == '-sParms':
        sParms = str(sys.argv[i+1]).split(',')
    if sys.argv[i] == '-fastSim':
        fastSim = True
    if sys.argv[i] == '-MCatNLO':
        MCatNLO = True
    if sys.argv[i] == '-isData':
        isData = True
    if sys.argv[i] == '-s':
        sample = sys.argv[i+1]

makeCrabConfig()

#if os.path.exists(cmsswSkelFile) == False:
#    print 'CMSSW skeleton file does not exist. Exiting'
#    sys.exit()
#
#
##print '\nGetting global tag from DBS...'
#if isData == True :
#    print 'Running on data sample.'
#
#if( global_tag_flag != '' ):
#	print '\nUsing \'' + global_tag_flag + '\' specified by -gtag flag.\n'
#	global_tag = global_tag_flag
#        if sParms > 0:
#            print 'Including sParmMaker with parameters %s.\n'%sParms
#	makeCMSSWConfig(cmsswSkelFile)
#	makeCrabConfig()
#else :
#    global_tag = '';
#    dbs_result = '';
#    if dbs_url == '' :
#      command = 'dbsql find config.name,config.content where dataset=' + dataSet.replace("AODSIM", "GEN-SIM-RECO") + '>config.content; while read line; do globaltag=`echo $line | sed -n \'s/^.*process.GlobalTag.globaltag = \([^p]*\).*$/\\1/p\'`; if [ "$globaltag" != "" ]; then echo $globaltag; break; fi; done <config.content; rm config.content';
#    else:
#      command = 'python $DBSCMD_HOME/dbsCommandLine.py -c search --url=' + dbs_url + ' --query="find config.name,config.content where dataset=' + dataSet.replace("AODSIM", "GEN-SIM-RECO") + '">config.content; while read line; do globaltag=`echo $line | sed -n \'s/^.*process.GlobalTag.globaltag = \([^p]*\).*$/\\1/p\'`; if [ "$globaltag" != "" ]; then echo $globaltag; break; fi; done <config.content; rm config.content';
#    #print command
#
#    length = len( os.popen(command).readlines() )
#    if( length > 0 ):
#      lines = os.popen(command);
#      for i in lines.readlines():
#        dbs_result = re.sub('\n', '', i)
#        global_tag = re.sub('#.+$', '', dbs_result)
#        if( global_tag != '' and global_tag_flag == '' ):
#            print '\nDBS Query results:\t\'' + dbs_result + '\' ?\n'
#            print 'Use global tag from DBS:\t\'' + global_tag + '\' ?\n'
#            answer = raw_input('[y/n]?')
#            while(answer != 'y' and answer != 'n'): 
#                print 'Please pick either \'y\' or \'n\''
#                answer = raw_input('[y/n]?')
#            if answer == 'n':
#                print 'Enter alternative Global Tag:'
#                global_tag = raw_input('new global tag:')
#            if sParms > 0:
#                print 'Including sParmMaker with parameters %s.\n'%sParms
#            makeCMSSWConfig(cmsswSkelFile)
#            makeCrabConfig()
#    else: 
#      print '\nGlobal tag not found in DBS. Use -gtag to set global tag. Exiting...\n'
#      sys.exit()
