#! /usr/bin/env python
# Sept. 26, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: sailfish, picard
# Uses ~ 5GB RAM


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


## TODO Hardcoded to starcluter-cloudbiolinux - need to generalize eventually using config file
workerUserName = 'ubuntu'
sailfishPath = '/home/ubuntu/bin/Sailfish-0.6.3-Linux_x86-64'
picard = '/usr/local/share/java/picard-1.96/'
headNFSPath = '/home/cbl_data/'
localWDPath = '/mnt/' # Execution node wd on cloudbiolinux
subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/'+workerUserName+'/',headPath=headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'QUANTsailfish'
wd = os.path.join(localWDPath, evalName)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))
syn = synapseclient.Synapse()
syn.login()


## Get index
#index = sr.copyRefToWorkerNode(args.idx,headNFSPath,wd, syn)
index = os.path.join(headNFSPath, 'Hsapiens_Gencode19_sailfish')


## Get submission
submission = syn.getSubmission(args.bam, downloadFile = False)
localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localBAMfilePath		



## Run samtofastq step
R1file = os.path.join(wd, prefix + '_R1.fastq')
R2file = os.path.join(wd, prefix + '_R2.fastq')
outputFile = os.path.join(wd, prefix + '.samtofastq')

cmd = ' '.join(['java -Xmx2G -jar', os.path.join(picard, 'SortSam.jar'), 'INPUT=', localBAMfilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(picard, 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file, 'TMP_DIR=', os.path.join(wd, 'tmp')])

if not os.path.exists(R1file) and not os.path.exists(R2file):
	print '%s' % time.asctime()
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print '%s' % time.asctime()
	print '\n\nFinished SamToFastq\n'



## Run SAILFISH step
#os.environ["PATH"] = ':'.join([os.path.join(sailfishPath, 'bin'), '$PATH'])
#os.environ["LD_LIBRARY_PATH"] = ':'.join([os.path.join(sailfishPath, 'lib'), '$LD_LIBRARY_PATH'])

outputDir = os.path.join(wd, prefix+'_quant')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

cmd = ' '.join(['sailfish quant -i', index, '-l "T=PE:O=><:S=SA" -1', R1file, '-2', R2file, '-o', outputDir])
if not os.path.exists(os.path.join(wd, prefix+'_quant.sf')) and not os.path.exists(os.path.join(outputDir, 'quant.sf')):
	print '%s' % time.asctime()
	print >> commandsFile, '%s' % cmd
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
act = Activity(name='transcript quantitation', description='Alignment-free transcript quantitation using Sailfish.', executed=['syn2325155', cf.id])
act.used(target=submission.entityId, targetVersion=submission.versionNumber)
act.used(args.idx)

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




# clean up BAM files and FASTQ
os.remove(localBAMfilePath)
os.remove(R1file, R2file)

# change status of BAM 
# status = syn.getSubmissionStatus(submission)
# status.status = 'SCORED' # Scored is functioning as "finished" for now.
# status = syn.store(status)
