#! /usr/bin/env python
'''
# Run generic_loader -h for usage information.
# KKD for Sage Bionetworks
# Aug. 16, 2013
'''

import argparse, synapseclient, os, sys
from synapseclient import Project, Folder, File, Activity, Evaluation

import seq_loading as sl

parser = argparse.ArgumentParser(description='Example generic bulk loader for seq workflows. Contains specific functionality for TopHat that can be commented out for other aligners.')
parser.add_argument('--params', dest='params', required=True, help='Config file.')
parser.add_argument('--password', dest='pwd', required=True, help='Your login password for Synapse')
args = parser.parse_args()

allArgDict = sl.readConfig(configFile=args.params)

#  Log in to Synapse and set up project
syn = synapseclient.Synapse()
syn.login(allArgDict['login'], args.pwd)
project = sl.setUpSynapseProject(allArgDict['folderlist'], pid=allArgDict['projectID'], syn=syn)


# Populate the synapse project
file_list = os.listdir(allArgDict['dataPath']) 
	for file in file_list:	
		if file.endswith('bam'):
			stepIDs = dict()

			for item in range(1,finalStep):
				step = str(item)
				filename = file.rstrip('.bam') + allArgDict[step].get('suffix')
				filePath = os.path.join(allArgDict[step].get('path'), filename)
				(th_version, specific_command) = sl.parseTophatCall(filePath)
				dataFolderID = existingFolders.get(allArgDict[step].get('folder'))
				id = sl.add_workflow_step_to_synapse(filePath,step=step,parentid=dataFolderID,cmd=specific_command, software=th_version, stepDict=allArgDict[step], syn=syn)	
				stepIDs[step] = id
