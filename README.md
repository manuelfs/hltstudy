hltstudy
==========

Code to make flat ntuples with HLT quantities.
It's been tested on `CMSSW_7_2_0_pre6`. 

#### Making ntuples
Issue the following commands on a machine with SLC6 and CMSSW_7_2_0_pre6:

    cmsrel CMSSW_7_2_0_pre6
    cd CMSSW_7_2_0_pre6/src
    cmsenv
    git cms-merge-topic --unsafe cms-l1t-offline:Stage1plusHLT_V0 ### Only for 2015 L1 menu
    git cms-addpkg L1TriggerConfig/L1GtConfigProducers ### Only for 2015 L1 menu
    cp -a /afs/cern.ch/user/g/ghete/public/L1Menu/L1Menu_Collisions2015_25ns_v1/xml/Rev.1.1/L1Menu_Collisions2015_25ns_v1_L1T_Scales_20101224_Imp0_0x102f.xml L1TriggerConfig/L1GtConfigProducers/data/Luminosity/startup/ ### Only for 2015 L1 menu
    git clone git@github.com:manuelfs/hltstudy
    scram b -j 4
    cmsRun hltstudy/babymaker/python/run_el15.py

This will create a flat ntuple named ntuple_hlt_el15.root in the
current directory with events that pass a trigger with HT>200 GeV
and pT of electron > 15 GeV.

In case of not having configured SSH in git, you can also check out the 
CfANtupler package with the http protocol

    git clone http://github.com/manuelfs/hltstudy

#### HLT producers
The file `setup_cff.py` was obtained with
    edmConfigFromDB --cff --configName /dev/CMSSW_7_1_1/AlternativeTrackingScenarios/GRun_TK1B/V28 --nopaths --services -PrescaleService > setup_cff.py

The files `hltstudy/babymaker/python/run_xxxx.py` are based in confDB paths that were obtained with

    hltGetConfiguration /users/manuelf/RA4ucsb/V5 --full --offline --mc --unprescale --process reHLT --globaltag auto:startup_GRun > run_hlt.py

The files `hltstudy/babymaker/python/cfg711/run_hlt_confdb_x.py` were obtained with

    hltGetConfiguration /users/jaehyeok/HLTRun2/RA4/V7 --full --offline --mc --unprescale --process reHLT --globaltag auto:startup_GRun > run_hlt.py
