hltstudy
==========

Code to make flat ntuples with HLT quantities.
It's been tested on `CMSSW_7_4_0_pre9`. 

#### Making ntuples
Issue the following commands on a machine with SLC6 and CMSSW_7_4_0_pre9:

    cmsrel CMSSW_7_4_4
    cd CMSSW_7_4_4/src
    cmsenv
    git cms-addpkg HLTrigger/Configuration
    git cms-checkdeps -A -a
    git clone git@github.com:manuelfs/hltstudy
    scram b -j$(getconf _NPROCESSORS_ONLN)

    cd hltstudy/babymaker/python/
    ./get_config.py /users/manuelf/CMSSW_7_4_X/ForBabies --output cfg_el15_740pre9.py
    cmsRun cfg_el15_740pre9.py

This will create a flat ntuple named ntuple_hlt_el.root in the
current directory with events that pass a trigger with HT>200 GeV
and pT of electron > 15 GeV.

In case of not having configured SSH in git, you can also check out the 
CfANtupler package with the http protocol

    git clone http://github.com/manuelfs/hltstudy

#### HLT producers
The file `setup_cff.py` was obtained with

    edmConfigFromDB --cff --configName /dev/CMSSW_7_4_0/GRun --nopaths --services -PrescaleService,-EvFDaqDirector,-FastMonitoringService > setup_cff.py

The files `hltstudy/babymaker/python/cfg_xxxx.py` are based in confDB paths that were obtained with

    ./get_config.py /users/manuelf/CMSSW_7_4_X/ForBabies --output cfg_el15_740.py
