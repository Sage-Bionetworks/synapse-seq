#! /usr/bin/env python
# Sept. 26, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: sailfish, picard


import os, argparse, subprocess, sys, time
import synapseclient, synapseseq
from synapseclient import Activity, File
from synapseseq import seq_running as sr

parser = argparse.ArgumentParser(description='Quantifies transcripts in input BAM using sailfish.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--index', dest='idx', required=True, help='Synapse ID of transcript index required by sailfish.', default=None)
parser.add_argument('--output', dest='out', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


## TODO Hardcoded to cloudbiolinux - need to generalize eventually using config file
sailfishPath = ''
picard = ''
headNFSPath = '/mnt/transient_nfs/'
localWDPath = '/mnt/galaxyData' # Execution node wd on cloudbiolinux
subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/ubuntu/',headPath=headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'QUANTsailfish'
wd = os.path.join(localWDPath, evalName)
if not os.path.exists(wd):
	os.mkdir(wd)
syn = synapseclient.Synapse()
syn.login()


## Get index
index = copyRefToWorkerNode(args.idx,headNFSPath,wd)


## Get submission
submission = syn.getSubmission(args.bam, downloadFile = False)
localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localBAMfilePath		



## Run samtofastq step
R1file = os.path.join(prefix + '_R1.fastq')
R2file = os.path.join(prefix + '_R2.fastq')
outputFile = os.path.join(prefix + '.samtofastq')
cmd = ' '.join(['java -Xmx2G -jar', os.path.join(picard, SamToFastq.jar), 'INPUT=', localBAMfilePath, 'FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file])
if not os.path.exists(R1file) and not os.path.exists(R2file):
	print '%s' % time.asctime()
	print >> comandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print '%s' % time.asctime()
	print '\n\nFinished SamToFastq\n'


## Run SAILFISH step
os.environ["PATH"] = ':'.join([os.path.join(sailfishPath, 'bin'), '$PATH'])
os.environ["LD_LIBRARY_PATH"] = ':'.join([os.path.join(sailfishPath, 'lib'), '$LD_LIBRARY_PATH'])

outputDir = os.path.join(wd, prefix+'_quant')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

cmd = ' '.join(['sailfish quant -i', index, '--reads', R1file, R2file, '-o', outputDir, "-l T=PE:O=<>:S=U"])
if not os.path.exists(os.path.join(wd, prefix+'_quant.sf')) and not os.path.exists(os.path.join(outputDir, 'quant.sf')):
	print '%s' % time.asctime()
	print >> comandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print '%s' % time.asctime()
	
	os.rename(os.path.join(outputDir, 'quant.sf'), os.path.join(wd, prefix+'_quant.sf'))
	os.rename(os.path.join(outputDir, 'quant_bias_corrected.sf'), os.path.join(wd, prefix+'_quant_bias_corrected.sf'))
	os.rename(os.path.join(outputDir, 'reads.count_info'), os.path.join(wd, prefix+'_count_info.sf'))



## Load results to synapse

# Set up provenance.
print 'Loading %s to Synapse.' % cfPath
commandsFile.close()
cf = File(path=cfPath, description='Job commands.', parentId=args.out, synapseStore=True)
cf = syn.store(cf, activityName='quant_evaluation', executed=['https://github.com/Sage-Bionetworks/synapse-seq/blob/master/scripts/eval_quant_sailfish.py'])
act = Activity(name='transcript quantitation', description='Alignment-free transcript quantitation using Sailfish.', executed=['syn2325155', cf.id], used=[args.idx, {'reference':{'target':submission.entityId, 'targetVersion':submission.versionNumber}, 'wasExecuted':False}])

# Load raw quant file
print 'Loading %s to Synapse.' % os.path.join(wd, prefix+'_quant.sf')
quantEntity = File(path=os.path.join(wd, prefix+'_quant.sf'), name=prefix+'_quant.sf', description='Quantified transcript isoforms.', parentId=args.out, synapseStore=True)	
quantEntity = syn.store(quantEntity, forceVersion=False, activity=act)
syn.setAnnotations(quantEntity, annotations=dict(fileType='quantitation',normalized='TPM',summaryLevel='transcript',biasCorrection='False'))
print 'new entity id %s' % quantEntity.id

# Load bias-corrected file.
print 'Loading %s to Synapse.' % os.path.join(wd, prefix+'_quant_bias_corrected.sf')
quantBCEntity = File(path=os.path.join(wd, prefix+'_quant_bias_corrected.sf'), name=prefix+'_quant_bias_corrected.sf', description='Quantified transcript isoforms.', parentId=args.out, synapseStore=True)	
quantBCEntity = syn.store(quantBCEntity, forceVersion=False, activity=act)
syn.setAnnotations(quantBCEntity, annotations=dict(fileType='quantitation',normalized='TPM',summaryLevel='transcript',biasCorrection='True'))
print 'new entity id %s' % quantBCEntity.id




# clean up BAM files
os.remove(localBAMfilePath)

# change status of BAM 
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)
