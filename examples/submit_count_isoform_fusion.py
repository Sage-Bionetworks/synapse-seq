#! /usr/bin/env python
# Oct. 7, 2013
# KKD for Sage Bionetworks
# Submit data files to an evaluation

import synapseclient, os, argparse, os.path
from synapseclient import Evaluation, Submission, File, Entity, Folder, Project
import seq_loading as sl

parser = argparse.ArgumentParser(description='Runs XYZ pipeline via synapse.')
parser.add_argument('--bamdir', dest='bam', required=True, help='Location of BAM files to process.')
parser.add_argument('--projectID', dest='pid', required=False, help='Synapse ID of project to which files should be loaded in Synapse -- if does not exist, then specify project name instead.', default=None)
parser.add_argument('--projectName', dest='pname', required=False, help='Project name to which files should be loaded in Synapse.', default=None)
args = parser.parse_args()

syn = synapseclient.Synapse()
syn.login()


(project, foldersDict) = sl.setUpSynapseProject(list(['BAM', 'counts', 'fusions']), syn, pid=args.pid, pname=args.pname)


countEval = syn.getEvaluation(2275607)
isoformEval = syn.getEvaluation(2275608)
fusionEval = syn.getEvaluation(2275609)
#tmp = syn.joinEvaluation(countEval)
#tmp = syn.joinEvaluation(isoformEval)
#tmp = syn.joinEvaluation(fusionEval)
profile = syn.getUserProfile()

# Record BAM files in synapse
file_list = os.listdir(args.bam) 
for file in file_list:	
	if file.endswith('bam'):
		filePath = os.path.join(args.bam, file)
#		print '%s' % filePath
		newFile = File(filePath, description = 'BAM file of aligned reads.', parentId = foldersDict['BAM'], synapseStore = False)
		newFile = syn.store(newFile)
	
		# Submit all BAM files to count evaluation, isoform eval, and fusion eval
		submission = syn.submit(entity=newFile, evaluation = countEval,  name = 'submissionTest', teamName = profile['displayName'])
		print 'Submitted %s to %s' % (newFile.name, countEval.name)
		submission = syn.submit(entity=newFile, evaluation = isoformEval,  name = 'submissionTest', teamName = profile['displayName'])
		print 'Submitted %s to %s' % (newFile.name, isoformEval.name)
		submission = syn.submit(entity=newFile, evaluation = fusionEval,  name = 'submissionTest', teamName = profile['displayName'])
		print 'Submitted %s to %s' % (newFile.name, fusionEval.name)
#submission = syn.submit(entity=foldersDict['BAM'], evaluation = syn.getEvaluation(evalID),  name = 'submissionTest', teamName = profile['displayName'])


### Notes:
# Can I submit a folder to an evaluation? 
# -->> NO....
# How are evaluations exposed to users? i.e. how will they know the eval id? 
# -->> They won't, but they might know the name.
# How to break up data between projects? i.e. could data be generated for more than one project within the window that the cron job runs?