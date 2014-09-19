#! /usr/bin/env python

import subprocess
import os

def time(input_name, prefix, path_name):
    output_name = os.path.splitext(os.path.basename(file_name))[0]
    subprocess.call(['./get_config.py',
                     '--paths',
                     'all',
                     '--max-events',
                     '-1',
                     '--timing',
                     '--input',
                     input_name,
                     path_name])
    subprocess.call(['cmsRun',
                     'config.py'])
    subprocess.call(['mv',
                     'DQM_V0001_R000000001__HLT__FastTimerService__All.root',
                     prefix+'_'+output_name+'.root'])

file_names = [
    'file:/nfs-7/userdata/ald77/HLT/qcd_300_470.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_1400_1800.root',
    'file:/nfs-7/userdata/jaehyeok/HLT/ttbar-13tev-850evts.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_5_10.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_10_20.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_20_30.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_30_50.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_50_80_.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_80_120.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_120_170.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_170_300.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_470_600.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_600_800.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_800_1000.root',
    'file:/nfs-7/userdata/ald77/HLT/qcd_1000_1400.root',
    ]

#file_names = [
#    'file:/nfs-7/userdata/jaehyeok/HLT/ttbar-13tev-850evts.root',
#    ]

for file_name in file_names:
    time(file_name, 'jae', '/users/jaehyeok/HLTRun2/SingleLeptonSUSY')
    time(file_name, 'manuel', '/users/manuelf/SingleLeptonSUSY')
    time(file_name, 'olivito', '/users/olivito/dev_7_1_0/htmet')
