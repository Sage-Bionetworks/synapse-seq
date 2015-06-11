#! /usr/bin/env python
# June 2, 2015
# Kristen K Dang for Sage Bionetworks
# Requires: rnastar, picard
# Uses ~ xGB RAM


import os, argparse, subprocess, sys, time, shutil
import synapseclient, synapseseq
from synapseclient import Activity, File
#from synapseclient.table import Table
from synapseseq import seq_running as sr
from synapseseq import synseqConfig as syncfg


parser = argparse.ArgumentParser(description='Aligns reads in input BAM using rna-star.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--ref', dest='ref', required=True, help='Synapse ID of reference genome used for alignment.', default=None)
parser.add_argument('--output', dest='out', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--threads', dest='threads', required=True, help='Number cpus to use per job')
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', syncfg.localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/'+syncfg.workerUserName+'/',headPath=syncfg.headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'ALIGNstar'
wd = os.path.join(syncfg.localWDPath, evalName)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))
syn = synapseclient.Synapse()
syn.login()


## Get index
#index = sr.copyRefToWorkerNode(args.idx,syncfg.headNFSPath,wd, syn)
index = os.path.join(syncfg.headNFSPath, syncfg.GRCh38_star_RL100)


## Get submission
submission = syn.getSubmission(args.bam, downloadFile = False)
BAMentity = syn.get(entity=submission.entityId, version=submission.versionNumber, downloadFile = True, downloadLocation = wd)
localBAMfilePath = os.path.join(wd, submission.name)
#localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localBAMfilePath		



## Run samtofastq step
## ** Need to handle SE reads
R1file = os.path.join(wd, prefix + '_R1.fastq')
R2file = os.path.join(wd, prefix + '_R2.fastq')
outputFile = os.path.join(wd, prefix + '.samtofastq')


cmd = ' '.join(['java -Xmx2G -jar', os.path.join(syncfg.picard, 'SortSam.jar'), 'INPUT=', localBAMfilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(syncfg.picard, 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file, 'TMP_DIR=', os.path.join(wd, 'tmp'), 'VALIDATION_STRINGENCY=SILENT'])

if not os.path.exists(R1file) and not os.path.exists(R2file):
	print 'samtofastq start %s' % time.asctime()
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print 'samtofastq end %s' % time.asctime()



## Run star step
outputDir = os.path.join(wd, prefix+'_align')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

cmd = ' '.join(['STAR --runMode alignReads --runThreadN', args.threads, '--genomeDir', index, '--readFilesIn', R1file, R2file, '--outFileNamePrefix', os.path.join(outputDir,prefix+'.'), '--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within'])
if not os.path.exists(os.path.join(outputDir, prefix+'.Aligned.sortedByCoord.out.bam')):
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	

## Load results to synapse

# Set up provenance.
print 'Loading %s to Synapse.' % cfPath
commandsFile.close()
cf = File(path=cfPath, description='Job commands.', parentId=args.out, synapseStore=True)
cf = syn.store(cf, activityName='alignment_evaluation', executed=['https://github.com/Sage-Bionetworks/synapse-seq/blob/master/scripts/eval_align_star.py'])
act = Activity(name='align to human', description='Spliced alignment of RNA reads to human genome using Star.', executed=['syn2243150', cf.id]) 
act.used(target=submission.entityId, targetVersion=submission.versionNumber)
act.used(args.ref)

# Load BAM file
print 'Loading %s to Synapse.' % os.path.join(outputDir,prefix+'.Aligned.sortedByCoord.out.bam')
alignEntity = File(path=os.path.join(outputDir,prefix+'.Aligned.sortedByCoord.out.bam'), name=prefix+'.Aligned.sortedByCoord.out.bam', description='Aligned reads in BAM format.', parentId=args.out, synapseStore=True)	
alignEntity = syn.store(alignEntity, forceVersion=False, activity=act)
syn.setAnnotations(alignEntity, annotations=dict(fileType='BAM',includeUnmapped=True))
print 'new entity id %s' % alignEntity.id
 



## Use this code after related JIRA is resolved
## Load metrics to Synapse table
# table = syn.get(syncfg.sailfishMetricsTable)
# metrics = list()
# metrics.append(prefix)
# with open(os.path.join(wd, prefix+'_count_info.sf')) as metricsFile:	
# 	for line in metricsFile:
# 		metrics.append(line.split()[1])
# metricsFile.closed
# syn.store(Table(table, [metrics]))



# clean up local files: BAM, FASTQ, sailfish files
os.remove(localBAMfilePath)
os.remove(R1file)
os.remove(R2file)
shutil.rmtree(wd)

# change status of BAM 
# status = syn.getSubmissionStatus(submission)
# status.status = 'SCORED' # Scored is functioning as "finished" for now.
# status = syn.store(status)
