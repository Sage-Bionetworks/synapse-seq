#! /usr/bin/env python
# Nov 24, 2015
# Kristen K Dang for Sage Bionetworks
# Requires: featurecounts



import os, argparse, subprocess, sys, shutil, yaml
import synapseclient, synapseseq
from synapseclient import File, Table, Activity
from synapseseq import seq_running as sr



parser = argparse.ArgumentParser(description='Submits a workflow job to a grid queue manager such as SGE.')
parser.add_argument('submission', help='Submission id of input file.')
parser.add_argument('--eval', required=True, help='Synapse ID of the eval to run.')
parser.add_argument('--conf', required=True, help='Path to local config file.')
parser.add_argument('--debug', required=False, help='Print all commands to output file even if not executed.', action='store_true')
args = parser.parse_args()

## Load workflow params
config = yaml.load(file(args.conf))


## Set up worker node
subprocess.call(' '.join(['sudo chown', config['system']['workerUserName']+':'+config['system']['workerUserName'], config['system']['localWDPath']]), shell = True)


## Log in to Synapse
sr.getSynConfigFromHead(localPath='/home/'+config['system']['workerUserName']+'/',headPath=config['system']['headNFSPath'])
syn = synapseclient.login()
evaluation = syn.getEvaluation(args.eval)


## Set up output directory
wd = os.path.join(config['system']['localWDPath'], evaluation.name)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))



## Get submission
submission = syn.getSubmission(args.submission, downloadFile = False)
BAMannotations = syn.getAnnotations(entity=submission.entityId)
localInputFilePath = os.path.join(wd, submission.name)
if not os.path.exists(localInputFilePath):
	localInputFilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=BAMannotations['bucket'],extKey=BAMannotations['keyname'])
# 	submission = syn.getSubmission(args.submission, downloadFile = True, downloadLocation = wd)
# 	if not localInputFilePath == submission.filePath:
# 		localInputFilePath = submission.filePath
prefix = os.path.basename(localInputFilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evaluation.name, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localInputFilePath		



submission = syn.getSubmission(args.bam, downloadFile = False)
localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print '%s' % localBAMfilePath	



## Get GTF
gtfEntity = syn.get(config['featurecounts']['gtf'], downloadLocation = config['system']['localWDPath'])
if gtfEntity.name.endswith('.gz'):
	subprocess.call('gunzip '+gtfEntity.path, shell = True)
	gtfName = gtfEntity.name.strip('.gz')
else:
	gtfName = gtfEntity.name
gtf = os.path.join(config['system']['localWDPath'], gtfName)


### Run featurecounts
outCountsFile = os.path.join(wd, '_'.join([prefix, config['featurecounts']['count-unit'], '_counts.txt']))
os.chdir(wd) # Attempt to redirect temp output files to this location, since they seem to go to the wd with no option to specify elsewhere.


cmd = ' '.join(['featureCounts', config['featurecounts']['commandParams'], '-a', gtf, '-o', outCountsFile, outBAMfile])

if not os.path.exists(outCountsFile):
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
elif args.debug is True:
	print >> commandsFile, '%s' % cmd

## featurecounts
table = syn.get(config['featurecounts']['metricsTable'])
metrics = [prefix, args.submission]
outLogFile = os.path.join(wd, '_'.join([prefix, config['featurecounts']['count-unit'], '_counts.txt.summary']))
with open(outLogFile) as metricsFile:	
	for line in metricsFile:
		if line.startswith('Status'): continue
		metrics.append(line.split()[1])
metricsFile.close()
syn.store(Table(table, [metrics]))


	
commandsFile.close()


# Load results to synapse

## Set up provenance.
print 'Loading %s to Synapse.' % cfPath
cf = File(path=cfPath, description='Job commands.', parentId=config['workflow']['output'], synapseStore=True)
cf = syn.store(cf)

## Load output counts file 
print 'Loading %s to Synapse.' % outCountsFile
outCountsEntity = File(path=outCountsFile, name=os.path.basename(outCountsFile), description='Gene or exon counts based on alignment to human genome.', parentId=config['featurecounts']['output'], synapseStore=True)	

outCountsEntity = syn.store(outCountsEntity, forceVersion=False, activityName='Quantitation', activityDescript='Counting reads that align to gene models.', used=[outBAMEntity, config['featurecounts']['gtf']], executed=[cf.id, 'syn2807330'])

syn.setAnnotations(outCountsEntity, annotations=config['featurecounts']['annotations'])
print 'new entity id %s' % outCountsEntity.id

 

## clean up local files
if 'keep-local' not in config['workflow']:
	os.remove(localInputFilePath)
	shutil.rmtree(wd)

## change status of input submission
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)

