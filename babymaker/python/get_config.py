#! /usr/bin/env python

import argparse
import subprocess

def WriteConfig(file):
    file.write("process.load('hltstudy.babymaker.setup_cff')\n")
    file.write("process.load('hltstudy.babymaker.babymakermu_cfi')\n")
    file.write("\n")
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

def WriteOutput(file):
    file.write("# Output\n")
    file.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
    file.write("process.out = cms.OutputModule(\n")
    file.write("        'PoolOutputModule',\n")
    file.write("        fileName     = cms.untracked.string('ntuple_hlt_mu15.root'), ##\n")
    file.write("        dropMetaData = cms.untracked.string('ALL')\n")
    file.write("        #SelectEvents = cms.untracked.PSet(\n")
    file.write("            #SelectEvents = cms.vstring('HLT_L1HTT')\n")
    file.write("        #),\n")
    file.write("\n")
    file.write(")\n")
    file.write("process.outpath      = cms.EndPath(process.out)\n")
    file.write("\n")
    file.write("process.out.outputCommands = cms.untracked.vstring( 'drop *' ) ##\n")
    file.write("process.out.outputCommands.extend(cms.untracked.vstring('keep *_*babymaker*_*_*')) ##\n")
    file.write("\n")
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

def WriteSchedule(file, path_list):
    seq_lines = []
    path_lines = []
    names = []

    for path in path_list:
        path_name = GetPathName(path)
        path_def = GetPathDef(path)

        seq_line = 'process.Seq' + path_name + ' = cms.Sequence( ' + path_def + ' )\n'
        path_line = 'process.' + path_name + ' = cms.Path( process.Seq' + path_name + ' + process.reco + process.babymaker)\n'

        names.extend([path_name])
        seq_lines.extend([seq_line])
        path_lines.extend([path_line])

    for seq in seq_lines:
        file.write(seq)
    file.write("\n")

    for path in path_lines:
        file.write(path)
    file.write("\n")

    if len(names) > 0:
        file.write('process.schedule = cms.Schedule( process.' + names[0] + ', process.outpath ) ##\n')            
        file.write("\n")

def WriteCustomization(file):
    file.write("##\n")
    file.write("from SLHCUpgradeSimulations.Configuration.postLS1Customs import *\n")
    file.write("process = customise_HLT( process )\n")

parser = argparse.ArgumentParser(
    description = 'Gets a configuration file and adds RECO products to the output.'
    )

parser.add_argument('config',
                    help = 'A confDB HLT configuration path')
parser.add_argument('--output',
                    default = 'config.py',
                    help = 'File to with python configuration is written')
parser.add_argument('--input',
                    default = 'file:/nfs-7/userdata/jaehyeok/HLT/ttbar-13tev-850evts.root',
                    help = 'Path to input data file')
parser.add_argument('--max-events',
                    type = int,
                    default = -1,
                    dest = 'maxevents',
                    metavar = 'EVENTS',
                    help = 'Maximum number of events to run on (-1 for all events)')
                    

args = parser.parse_args()

base_file_name = subprocess.check_output('mktemp')[:-1]
base_file = open(base_file_name, 'w')
subprocess.call(['hltGetConfiguration',
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
                 '--globaltag',
                 'auto:startup_GRun'],
                stdout=base_file)
base_file.close()

base_file = open(base_file_name, 'r')
lines = base_file.readlines()

output_file = open(args.output, 'w')

path_list = []

for line in lines:
    if line.find('process.HLTConfigVersion = cms.PSet(') != -1:
        WriteConfig(output_file)
#    Uncommented these lines to make Manuel's muon configuration from ConfDB work
#    elif line.find('hltESPAK4PFL1L2L3') != -1:
#        line = line.replace('hltESPAK4PFL1L2L3','hltESPAK4PFCorrection')
#    elif line.find('hltESPAK4PFNoPUL1L2L3') != -1:
#        line = line.replace('hltESPAK4PFNoPUL1L2L3','hltESPAK4PFCorrection')
#    elif line.find('hltESPAK4CaloL1L2L3') != -1:
#        line = line.replace('hltESPAK4CaloL1L2L3','hltESPAK4CaloCorrection')
    elif line.find("process.HLT") != -1 and line.find("cms.Path") != -1:
        path_list.extend([line])
        line = ''
    elif line.find("process.source") != -1:
        WriteGenJets(output_file)
    elif line.find("# enable the TrigReport and TimeReport") != -1:
        WriteMixing(output_file)
    elif line.find("wantSummary = cms.untracked.bool( True )") != -1:
        line = line.replace("True","False")
    elif line.find('# override the GlobalTag, connection string and pfnPrefix') != -1:
        WriteOutput(output_file)
        WriteSchedule(output_file, path_list)

    output_file.write(line)

WriteCustomization(output_file)

base_file.close()
output_file.close()
