#! /usr/bin/env python
# Oct. 8, 2013
# KKD for Sage Bionetworks
# Finds new BAM files and submits them to an evaluation workflow.

import synapseclient, os, argparse, sys, os.path
from synapseclient import File, Entity, Folder, Project, Activity
import seq_loading as sl

parser = argparse.ArgumentParser(description='Submits to XYZ pipeline via synapse.')
parser.add_argument('--jobs', dest='jobs', required=False, help='Maximum number of concurrent jobs to run in parallel sections.', default = 2)
parser.add_argument('--output', dest='pid', required=True, help='Synapse project or folder to which to add results.', default = 'syn1972151')
parser.add_argument('--projectName', dest='pname', required=False, help='Project name to which files should be loaded in Synapse.', default=None)
parser.add_argument('--printout', dest='printout', required=False, help='Print jobs to run but do not run them.', action='store_true')
args = parser.parse_args()

syn = synapseclient.Synapse()
syn.login()

evalID = '2275607'
profile = syn.getUserProfile() ## If this is run from cron, need to chose another way to get the user profile.
snpEval = syn.getEvaluation(evalID)
#tmp = syn.joinEvaluation(snpEval)

### Detection/setup of project/folder needs work:
# 		Need to handle either folder or project creation.
#		What to do if project name is already taken or invalid
#		Do sub-evals check for and create their specific sub-folders?

# Set up project/folder if it doesn't exist
if args.pname not None:
	(project, foldersDict) = sl.setUpSynapseProject(list(['BAM']), syn, pname=args.pname)
# Otherwise, just make sure the BAM folder is there and get its id.
else: 
	(project, foldersDict) = sl.setUpSynapseProject(list(['BAM']), syn, pid=args.pid)



# Get dict of existing submissions
existingBAMDict = dict()
for submission in syn.getSubmissions(evalID):
	BAMentity = syn.get(submission.entityId, downloadFile = False)
	existingBAMDict[BAMentity.path] = BAMentity


# Check for new BAM and submit
BAMDir = os.listdir(args.bam)
for dir in BAMDir:
	filesList = os.listdir(os.path.join(args.bam, dir)):
	for file in filesList:
		if file.endswith('.bam'):
			if file not in existingBAMDict:		
### For qsub, the evaluation code for this submission will have to contain the provenance calls.
				filePath = os.path.join(args.bam, dir, file)
				newFile = File(filePath, description = 'BAM file of aligned reads.', parentId = foldersDict['BAM'], synapseStore = False)
				
				## Try to extract this code to a function in seq_loading
				act = Activity(name='Alignment',  description='Align reads to genome.')
				#act.executed(target='tophatid', targetVersion=version)

				newFile = syn.store(newFile, activity = act)
				submission = syn.submit(entity=newFile, evaluation = evalID,  name = 'submissionTest', teamName = profile['displayName'])
				print 'Submitted %s to %s' % (newFile.name, countEval.name)

				### Could have multiple other submissions here: count, fusion, 


## Notes
# Can evaluations be used or executed entities in an activity?
# -->> Prefer to use/execute the code behind the eval?
# How to break up data between projects? i.e. could data be generated for more than one project within the window that the cron job runs?
# --> ?? Don't know yet...hope it doesn't come to parsing the sample sheet. If cron is hard-coded to check a specific directory per project, that will work. Not so much if data is dumped into same directory generically for all projects.
# How are evaluations exposed to users? i.e. how will they know the eval id? 
# -->> They won't, but they might know the name.