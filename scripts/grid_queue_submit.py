#! /usr/bin/env python
# June 2, 2015
# Kristen K Dang for Sage Bionetworks
# Requires:
# Uses ~ xGB RAM


import os, argparse, subprocess, sys, time, shutil
import synapseclient, synapseseq
from synapseclient import Activity, File
#from synapseclient.table import Table
from synapseseq import seq_running as sr


parser = argparse.ArgumentParser(description='Submits a workflow job to a grid queue manager such as SGE.')
parser.add_argument('--submission', required=True, dest='subm', help='Submission id of input file.')
parser.add_argument('--evalID', required=True, help='Synapse ID of the eval to run.')
parser.add_argument('--config', dest='conf', required=True, help='Path to local config file.')
args = parser.parse_args()

## Load workflow params
config = yaml.load(file(args.conf))


## Set up worker node
subprocess.call(' '.join(['sudo chown' config['workerUserName']+':'+config['workerUserName'], config['localWDPath']]), shell = True)


## Log in to Synapse
sr.getSynConfigFromHead(localPath='/home/'+config['workerUserName']+'/',headPath=config['headNFSPath'])
syn = synapseclient.Synapse()
syn.login()
evaluation = syn.getEvalutation(args.evalID)


## Set up output directory
wd = os.path.join(config['localWDPath'], evaluation.name)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))


## Get reference index path on head node
index = os.path.join(config['headNFSPath'], config[evaluation.name['index']])


## Get submission
submission = syn.getSubmission(args.subm, downloadFile = True, downloadLocation = wd)
localInputFilePath = os.path.join(wd, submission.name)
# put some code here for getting files using AWS CLI for legacy data
#localInputFilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localInputFilePath).rstrip('.'+config[evaluation.name['input-suffix']])
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localInputFilePath		



## Run evaluation
#*****	#
return: synid of executable, output file path, annotations dictionary

commandsFile.close()

## Load results to synapse (opt)
# ***** synid of executable, output entity description, annotations
if 'output' in config[evaluation.name]:
	# Set up provenance.
	print 'Loading %s to Synapse.' % cfPath
	cf = File(path=cfPath, description='Job commands.', parentId=config[evaluation.name['output']], synapseStore=True)
	cf = syn.store(cf)
	act = Activity(name=evaluation.name, description=evaluation.description, executed=[xxexecutablexxxx, cf.id]) 
	act.used(target=submission.entityId, targetVersion=submission.versionNumber)
	act.used(args.ref)

	# Load output file 
	print 'Loading %s to Synapse.' % outFilePath
	outputEntity = File(path=outFilePath, name=os.path.basename(outFilePath), description=xxxxx, parentId=config[evaluation.name['output']], synapseStore=True)	
	outputEntity = syn.store(outputEntity, forceVersion=False, activity=act)
	syn.setAnnotations(outputEntity, annotations=xxxx)
	print 'new entity id %s' % outputEntity.id
 

## Load metrics to Synapse table
# table = syn.get(config['sailfishMetricsTable'])
# metrics = list()
# metrics.append(prefix)
# with open(os.path.join(wd, prefix+'_count_info.sf')) as metricsFile:	
# 	for line in metricsFile:
# 		metrics.append(line.split()[1])
# metricsFile.closed
# syn.store(Table(table, [metrics]))


## clean up local files
if 'keep-local' not in config[evaluation.name]:
	os.remove(localInputFilePath)
	shutil.rmtree(wd)

# change status of input submission
# status = syn.getSubmissionStatus(submission)
# status.status = 'SCORED' # Scored is functioning as "finished" for now.
# status = syn.store(status)


## Run following evaluations, if any
workflowSteps = config['workflow'['order']]
for i in range(len(workflowSteps)-1):
	if workflowSteps[i] is not evaluation.name: continue
	elif i < (len(workflowSteps-1):
		nextEval = workflowSteps[i+1]
	
	

# Submit any following jobs and run
followingSubmission = syn.submit(entity=outputEntity, evaluation = config[nextEval['id']],  name = outputEntity.name, teamName = submission.teamName)
print 'Submitted %s to %s' % (outputEntity.name, config[nextEval['id']])

# Call this same script, with different parameters
cmd = ' '.join(['synseq_grid_queue_submit.py --submission', followingSubmission.id, '--evalID', config[nextEval['id']]], '--config', args.config)
print '%s' % cmd
subprocess.call(cmd, shell = True)
