#! /usr/bin/env python

import argparse
import subprocess
import os

def WriteConfig(file, reco, timing):
    file.write("process.load('hltstudy.babymaker.setup_cff')\n")
    if not timing:
        file.write("process.load('hltstudy.babymaker.IsoMuonProducer_cfi')\n")
        file.write("process.load('hltstudy.babymaker.IsoElectronProducer_cfi')\n")
        if reco:
            file.write("process.load('hltstudy.babymaker.babymakerREPLACEreco_cfi')\n")
        else:
            file.write("process.load('hltstudy.babymaker.babymakerREPLACE_cfi')\n")
    file.write("\n")
    if reco:
        file.write("# import of standard configurations\n")
        file.write("process.load('Configuration.StandardSequences.Services_cff')\n")
        file.write("process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')\n")
        file.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
        file.write("process.load('Configuration.EventContent.EventContent_cff')\n")
        file.write("process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')\n")
        file.write("process.load('Configuration.StandardSequences.GeometryRecoDB_cff')\n")
        file.write("#process.load('Configuration.StandardSequences.MagneticField_38T_cff') #Might want this, but probably okay to omit\n")
        file.write("process.load('Configuration.StandardSequences.RawToDigi_cff')\n")
        file.write("process.load('Configuration.StandardSequences.L1Reco_cff')\n")
        file.write("process.load('Configuration.StandardSequences.Reconstruction_cff')\n")
        file.write("process.load('Configuration.StandardSequences.EndOfProcess_cff')\n")
        file.write("#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') #Might want this, maybe dangerous to omit\n")
        file.write("\n")
        file.write("process.esp_SiStripGainESProducer = cms.ESPrefer('SiStripGainESProducer')\n")
        file.write("process.esp_SiStripQualityESProducer = cms.ESPrefer('SiStripQualityESProducer')\n")
        file.write("process.esp_EcalTrigTowerConstituentsMapBuilder = cms.ESPrefer('EcalTrigTowerConstituentsMapBuilder','hltESPEcalTrigTowerConstituentsMapBuilder')\n")
        file.write("process.esp_HcalTopologyIdealEP = cms.ESPrefer('HcalTopologyIdealEP')\n")
        file.write("process.esp_TrackerGeometricDetESModule = cms.ESPrefer('TrackerGeometricDetESModule')\n")
        file.write("process.esp_MuonDetLayerGeometryESProducer = cms.ESPrefer('MuonDetLayerGeometryESProducer','hltESPMuonDetLayerGeometryESProducer')\n")
        file.write("process.esp_TrackerDigiGeometryESModule = cms.ESPrefer('TrackerDigiGeometryESModule')\n")
        file.write("process.esp_TrackerRecoGeometryESProducer = cms.ESPrefer('TrackerRecoGeometryESProducer','hltESPTrackerRecoGeometryESProducer')\n")
        file.write("process.esp_GlobalTrackingGeometryESProducer = cms.ESPrefer('GlobalTrackingGeometryESProducer','hltESPGlobalTrackingGeometryESProducer')\n")
        file.write("\n")

def WriteGenJets(file):
    file.write("## Gen jets\n")
    file.write("process.load('RecoJets.Configuration.GenJetParticles_cff')\n")
    file.write("process.load('RecoJets.Configuration.RecoGenJets_cff')\n")
    file.write("process.antiktGenJets = cms.Sequence(\n")
    file.write("    process.genParticlesForJetsNoNu*\n")
    file.write("    process.ak4GenJetsNoNu\n")
    file.write("    )\n")
    file.write("\n")

def WriteMixing(file):
    file.write("# Pileup mixing\n")
    file.write("process.mix.input.nbPileupEvents.averageNumber = cms.double(40.000000)\n")
    file.write("process.mix.bunchspace = cms.int32(25)\n")
    file.write("process.mix.minBunch = cms.int32(-12)\n")
    file.write("process.mix.maxBunch = cms.int32(3)\n")
    file.write("process.mix.input.fileNames = cms.untracked.vstring(['/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/001CB469-A91E-E311-9BFE-0025907FD24A.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/009CB248-A81C-E311-ACD8-00259073E4F0.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/009F81D5-B21C-E311-966C-BCAEC50971D0.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/00B5BB8C-A91E-E311-816A-782BCB1F5E6B.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/00B8F676-BA1C-E311-BA87-0019B9CABFB6.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/00DD7446-B51D-E311-B714-001E6739CEB1.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/021E1B53-101D-E311-886F-00145EDD7569.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/022A782D-A51C-E311-9856-80000048FE80.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/026FE678-BA1C-E311-BEF5-00D0680BF90A.root', '/store/mc/Fall13/MinBias_TuneA2MB_13TeV-pythia8/GEN-SIM/POSTLS162_V1-v1/10000/02A10BDE-B21C-E311-AB59-00266CF327C0.root'])\n")
    file.write("from Configuration.AlCa.GlobalTag import GlobalTag\n")
    file.write("process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')\n")
    file.write("\n")

def WriteOutput(file, reco, timing):
    file.write("# Output\n")
    file.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
    file.write("process.out = cms.OutputModule(\n")
    file.write("        'PoolOutputModule',\n")
    file.write("        fileName     = cms.untracked.string('ntuple_hlt_REPLACE.root'), ##\n")
    file.write("        dropMetaData = cms.untracked.string('ALL')\n")
    file.write("        #SelectEvents = cms.untracked.PSet(\n")
    file.write("            #SelectEvents = cms.vstring('HLT_L1HTT')\n")
    file.write("        #),\n")
    file.write("\n")
    file.write(")\n")
    file.write("process.outpath      = cms.EndPath(process.out)\n")
    file.write("\n")
    file.write("process.out.outputCommands = cms.untracked.vstring( 'drop *' ) ##\n")
    if not timing:
        file.write("process.out.outputCommands.extend(cms.untracked.vstring('keep *_*babymaker*_*_*')) ##\n")
    file.write("\n")

    if reco:
        file.write("# Path and EndPath definitions for RECO\n")
        file.write("process.raw2digi_step = cms.Path(process.RawToDigi)\n")
        file.write("process.L1Reco_step = cms.Path(process.L1Reco)\n")
        file.write("process.reconstruction_step = cms.Path(process.reconstruction)\n")
        file.write("process.endjob_step = cms.EndPath(process.endOfProcess)\n")
        file.write("\n")
        file.write("process.reco = cms.Sequence( process.RawToDigi + process.L1Reco + process.reconstruction + process.endOfProcess )\n")
        file.write("\n")

def GetPathName(path):
        stripped = path.strip()
        start_of_name = stripped.find('.')
        end_of_name = stripped.find(' ')
        if start_of_name != -1 and end_of_name != -1:
            return stripped[start_of_name+1:end_of_name]
        else:
            raise Exception("GetPathName could not find path name in " + stripped)
            return ''

def GetPathDef(path):
    stripped = path.strip()
    start = stripped.find('cms.Path( ') + 10
    end = stripped.find(' )')
    if start!= -1 and end != -1:
        raw_path = stripped[start:end]
        return stripped[start:end].replace('process.HLTEndSequence','process.antiktGenJets + process.HLTEndSequence')
    else:
        raise Exception("GetPathDef could not find path definition in " + stripped)
        return ''

def WriteSchedule(file, path_list, path_opt, reco, timing):
    seq_lines = []
    path_lines = []
    names = []

    for path in path_list:
        path_name = GetPathName(path)
        path_def = GetPathDef(path)

        if timing and path_name == 'HLTriggerFirstPath':
            continue

        if timing and path_opt == 'all':
            if (path_name == 'Timing_HLT_Mu15_IsoVVVL_IterTrk02_v1' or
                path_name == 'HLT_Mu15_IterTrk02_PFHT300_v1' or
                path_name == 'Test_HLT_Mu15_v1' or
                path_name == 'HLT_Ele15_IsoVVVL_v1' or
                path_name == 'HLT_Ele15_PFHT300_v1'):
                continue

        maybe_electron = path_name.find('Ele') != -1
        maybe_muon = path_name.find('Mu') != -1

        baby_maker = ''
        iso_producer = ''
        if not timing:
            baby_maker = ' + process.babymaker'
            if maybe_electron and not maybe_muon:
                iso_producer = ' + process.IsoElectronProducer'
            elif maybe_muon and not maybe_electron:
                iso_producer = ' + process.IsoMuonProducer'

        reco_producer = ''
        if reco:
            reco_producer = ' + process.reco'

        seq_line = 'process.Seq' + path_name + ' = cms.Sequence( ' + path_def + ' )\n'
        path_line = 'process.' + path_name + ' = cms.Path( process.Seq' + path_name + iso_producer + reco_producer + baby_maker +' )\n'

        names.extend([path_name])
        seq_lines.extend([seq_line])
        path_lines.extend([path_line])

    if timing:
        file.write('process.SeqHLTriggerFirstPath = cms.Sequence( process.hltGetRaw + process.hltGetConditions + process.hltBoolFalse )\n')
    for seq in seq_lines:
        file.write(seq)
    file.write("\n")

    if timing:
        file.write('process.HLTriggerFirstPath = cms.Path( process.SeqHLTriggerFirstPath )\n')

    to_schedule = [1]
    if path_opt == 'all':
        to_schedule = [i for i in range(1,len(names)+1)]
    elif path_opt == 'select':
        print 'Available paths:'
        for i in range(0,len(names)):
            print '\t'+str(i+1)+': \t'+names[i]
        print ''

        to_schedule = raw_input('Enter comma separated list of paths to schedule: ').split(',')
        if to_schedule == ['']:
            to_schedule = []
        else:
            to_schedule = [int(i) for i in to_schedule]

    lep_type = 'mu'
    if len(names) > 0:
        already_added = set()
        for path in to_schedule:
            if 0 < path <= len(names):
                if not path in already_added:
                    if names[path-1].find('Ele') != -1:
                        lep_type = 'el'
                    file.write(path_lines[path-1])
                    already_added.add(path)
                else:
                    print 'Warning: WriteSchedule has already added path #'+str(path)+' ('+names[path-1]+') and will not add it again.'
            else:
                print 'Warning: WriteSchedule does not know about path #'+str(path)+' and will skip it.'

    return lep_type

def WriteTimingOutput(file):
    file.write("# FastTimerServiceClient\n")
    file.write("process.fastTimerServiceClient = cms.EDAnalyzer( 'FastTimerServiceClient',\n")
    file.write("                                                 dqmPath = cms.untracked.string( 'HLT/TimerService' )\n")
    file.write("                                                 )\n")
    file.write("\n")
    file.write("# DQM file saver\n")
    file.write("process.dqmFileSaver = cms.EDAnalyzer( 'DQMFileSaver',\n")
    file.write("                                       convention        = cms.untracked.string( 'Offline' ),\n")
    file.write("                                       workflow          = cms.untracked.string( '/HLT/FastTimerService/All' ),\n")
    file.write("                                       dirName           = cms.untracked.string( '.' ),\n")
    file.write("                                       saveByRun         = cms.untracked.int32(1),\n")
    file.write("                                       saveByLumiSection = cms.untracked.int32(-1),\n")
    file.write("                                       saveByEvent       = cms.untracked.int32(-1),\n")
    file.write("                                       saveByTime        = cms.untracked.int32(-1),\n")
    file.write("                                       saveByMinute      = cms.untracked.int32(-1),\n")
    file.write("                                       saveAtJobEnd      = cms.untracked.bool(False),\n")
    file.write("                                       forceRunNumber    = cms.untracked.int32(-1),\n")
    file.write("                                       )\n")
    file.write("\n")
    file.write("process.TimingOutput = cms.EndPath( process.fastTimerServiceClient + process.dqmFileSaver )\n")
    file.write('\n')

def WriteCustomization(file):
    file.write("##\n")
    file.write("from SLHCUpgradeSimulations.Configuration.postLS1Customs import *\n")
    file.write("process = customise_HLT( process )\n")

parser = argparse.ArgumentParser(
    description = 'Gets a configuration file and adds RECO products to the output.',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

parser.add_argument('config',
                    nargs='?',
                    default='/users/jaehyeok/HLTRun2/SingleLeptonSUSY_720pre6',
                    help = 'A confDB HLT configuration path')
parser.add_argument('--output',
                    default = 'config.py',
                    help = 'File to python configuration to')
parser.add_argument('--input',
                    default = 'file:/nfs-7/userdata/jaehyeok/HLT/ttbar-13tev-850evts.root',
                    help = 'Path to input data file')
parser.add_argument('--max-events',
                    type = int,
                    default = 30,
                    dest = 'maxevents',
                    metavar = 'EVENTS',
                    help = 'Maximum number of events to run on, -1 for all events')
parser.add_argument('--paths',
                    metavar = 'PATH_OPT',
                    choices=['select','first','all'],
                    default='select',
                    help = 'Choose to schedule manually selected paths, the first path, or all paths; options are %(choices)s')
parser.add_argument('--oldL1',
                    action = 'store_true',
                    help = 'Use 2012 L1 menu')

mode = parser.add_mutually_exclusive_group()
mode.add_argument('--reco',
                  action='store_true',
                  help = 'Enable RECO')
mode.add_argument('--timing',
                  action='store_true',
                  help = 'Enable timing')

args = parser.parse_args()

base_file_name = subprocess.check_output('mktemp')[:-1]
intermediate_file_name = subprocess.check_output('mktemp')[:-1]

base_file = open(base_file_name, 'w')
hlt_args=['hltGetConfiguration',
          args.config,
          '--full',
          '--offline',
          '--mc',
          '--unprescale',
          '--process',
          'reHLT',
          '--max-events',
          str(args.maxevents),
          '--input',
          args.input,
          '--globaltag']
if args.oldL1:
    hlt_args.extend(['auto:startup_GRun'])
else:
    hlt_args.extend(['auto:upgradePLS1',
                     '--l1-emulator',
                     'stage1,gt',
                     '--l1Xml',
                     'L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml'])
if args.timing:
    hlt_args.extend(['--timing'])

for argo in hlt_args:
    print argo,
print
subprocess.call(hlt_args, stdout=base_file)
base_file.close()

base_file = open(base_file_name, 'r')
intermediate_file = open(intermediate_file_name, 'w')

path_list = []
lep_type = 'mu'
found_timer_service = False

lines = base_file.readlines()
for line in lines:
    if line.find('process.HLTConfigVersion = cms.PSet(') != -1:
        WriteConfig(intermediate_file, args.reco, args.timing)
    elif line.find("process.HLT") != -1 and line.find("cms.Path") != -1:
        path_list.extend([line])
        line = ''
    elif line.find("process.source") != -1:
        WriteGenJets(intermediate_file)
    elif line.find("# enable the TrigReport and TimeReport") != -1:
        if args.reco:
            WriteMixing(intermediate_file)
    elif line.find("wantSummary = cms.untracked.bool( True )") != -1:
        line = line.replace("True","False")
    elif line.find('# override the GlobalTag, connection string and pfnPrefix') != -1:
        WriteOutput(intermediate_file, args.reco, args.timing)
        lep_type = WriteSchedule(intermediate_file, path_list, args.paths, args.reco, args.timing)
    elif line.find('process.FastTimerService.dqmTimeRange') != -1:
        line = line.replace('1000.','2000.')
    elif line.find('process.FastTimerService.dqmPathTimeRange') != -1:
        line = line.replace('100.', '2000.')
    elif line.find('process.FastTimerService.dqmModuleTimeRange') != -1:
        line = line.replace('40.', '800.')
    elif line.find('FastTimerServiceClient') != -1:
        found_timer_service = True
    elif line.find('CMSSW version specific customizations') != -1:
        if args.timing and not found_timer_service:
            WriteTimingOutput(intermediate_file)
    elif line.find('customiseGlobalTag(') != -1:
        line = "    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V3A',conditions='TrackerAlignmentExtendedError_2011Realistic_v1_mc,TrackerAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonDTAPEObjectsExtended_v0_mc,DTAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonCSCAPEObjectsExtended_v0_mc,CSCAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalSamplesCorrelation_mc,EcalSamplesCorrelationRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseShapes_mc,EcalPulseShapesRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseCovariances_mc,EcalPulseCovariancesRcd,frontier://FrontierProd/CMS_CONDITIONS')\n"
    elif line.find("process.HLTSchedule") != -1:
        line = '# '+line

    intermediate_file.write(line)

intermediate_file.write("if hasattr(process, 'hltCsc2DRecHits'):\n")
intermediate_file.write('    process.hltCsc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")\n')
intermediate_file.write('    process.hltCsc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")\n')
intermediate_file.write("if hasattr(process, 'cscReEmulTriggerPrimitiveDigis'):\n")
intermediate_file.write('    process.cscReEmulTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")\n')
intermediate_file.write('    process.cscReEmulTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")\n')

if args.timing:
    intermediate_file.write("process.DQMStore = cms.Service('DQMStore')\n")

WriteCustomization(intermediate_file)

intermediate_file.close()
base_file.close()

intermediate_file = open(intermediate_file_name, 'r')
output_file = open(args.output, 'w')
lines = intermediate_file.readlines()
for line in lines:
    output_file.write(line.replace('REPLACE',lep_type))
output_file.close()
intermediate_file.close()

os.remove(intermediate_file_name)
os.remove(base_file_name)
