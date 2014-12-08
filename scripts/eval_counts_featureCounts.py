#! /usr/bin/env python
# Nov. 4, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: featurecounts
# Uses ~ 60MB RAM


import os, argparse, subprocess, sys, time
import synapseclient, synapseseq
from synapseclient import Activity, File
from synapseseq import seq_running as sr

parser = argparse.ArgumentParser(description='Quantifies transcripts in input BAM using sailfish.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--gtf', dest='gtf', required=True, help='Synapse ID of GTF defining gene models to use for counting.', default=None)
parser.add_argument('--output', dest='out', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


## TODO Hardcoded to starcluter-cloudbiolinux - need to generalize eventually using config file
workerUserName = 'ubuntu'
headNFSPath = '/home/cbl_data/'
localWDPath = '/mnt/' # Execution node wd on cloudbiolinux
subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/'+workerUserName+'/',headPath=headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'COUNTfeaturecounts'
wd = os.path.join(localWDPath, evalName)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))
syn = synapseclient.Synapse()
syn.login()


## Get GTF
gtfEntity = syn.get(args.gtf, downloadFile = False)
if gtfEntity.name.endswith('.gz'):
	gtfName = gtfEntity.name.strip('.gz')
else:
	gtfName = gtfEntity.name
gtf = os.path.join(headNFSPath, gtfName)


## Get submission
submission = syn.getSubmission(args.bam, downloadFile = False)
localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localBAMfilePath		



## Run featurecounts step
outputFile = os.path.join(wd, prefix+'_gene_counts.txt')

os.chdir(wd) # Attempt to redirect temp output files to this location, since they seem to go to the wd with no option to specify elsewhere.
cmd = ' '.join(['featureCounts -p -t exon -g gene_id -a', gtf, '-o', outputFile, localBAMfilePath])
print '%s' % time.asctime()
print >> commandsFile, '%s' % cmd
subprocess.call(cmd, shell = True)
print '%s' % time.asctime()
	


## Load results to synapse

# Set up provenance.
print 'Loading %s to Synapse.' % cfPath
commandsFile.close()
cf = File(path=cfPath, description='Job commands.', parentId=args.out, synapseStore=True)
cf = syn.store(cf, activityName='count_evaluation', executed=['https://github.com/Sage-Bionetworks/synapse-seq/blob/master/scripts/eval_counts_featurecounts.py'])
act = Activity(name='Read counting', description='Counting aligned reads to GTF features using featurecounts.', executed=['syn2807330', cf.id])
act.used(target=submission.entityId, targetVersion=submission.versionNumber)
act.used(args.gtf) 

# Load raw count file
print 'Loading %s to Synapse.' % outputFile
quantEntity = File(path=outputFile, name=prefix+'_gene_counts.txt', description='Read counts summarized at gene level.', parentId=args.out, synapseStore=True)	
quantEntity = syn.store(quantEntity, forceVersion=False, activity=act)
syn.setAnnotations(quantEntity, annotations=dict(fileType='count',normalized='no',summaryLevel='gene',biasCorrection='False'))
print 'new entity id %s' % quantEntity.id



# clean up BAM files
os.remove(localBAMfilePath)
os.remove(commandsFile)

# change status of BAM 
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)
