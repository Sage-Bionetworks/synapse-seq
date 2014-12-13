#! /usr/bin/env python
# Sept. 26, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: sailfish, picard
# Uses ~ 5GB RAM


import os, argparse, subprocess, sys, time, shutil
import synapseclient, synapseseq
from synapseclient import Activity, File
#from synapseclient.table import Table
from synapseseq import seq_running as sr
from synapseseq import synseqConfig as syncfg


parser = argparse.ArgumentParser(description='Quantifies transcripts in input BAM using sailfish.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--index', dest='idx', required=True, help='Synapse ID of transcript index required by sailfish.', default=None)
parser.add_argument('--output', dest='out', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--stranded', dest='strand', required=False, help='Strandedness of the reads (see sailfish docs for explanation of options).', default='U')
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', syncfg.localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/'+syncfg.workerUserName+'/',headPath=syncfg.headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'QUANTsailfish'
wd = os.path.join(syncfg.localWDPath, evalName)
if not os.path.exists(wd):
	os.mkdir(wd)
if not os.path.exists(os.path.join(wd, 'tmp')):
	os.mkdir(os.path.join(wd, 'tmp'))
syn = synapseclient.Synapse()
syn.login()


## Get index
#index = sr.copyRefToWorkerNode(args.idx,syncfg.headNFSPath,wd, syn)
index = os.path.join(syncfg.headNFSPath, 'Hsapiens_Gencode19_sailfish')


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

cmd = ' '.join(['java -Xmx2G -jar', os.path.join(syncfg.picard, 'SortSam.jar'), 'INPUT=', localBAMfilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(syncfg.picard, 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file, 'TMP_DIR=', os.path.join(wd, 'tmp')])

if not os.path.exists(R1file) and not os.path.exists(R2file):
	print 'samtofastq start %s' % time.asctime()
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print 'samtofastq end %s' % time.asctime()



## Run SAILFISH step
os.environ["PATH"] = ':'.join([os.path.join(syncfg.sailfishPath, 'bin'), '$PATH'])
os.environ["LD_LIBRARY_PATH"] = ':'.join([os.path.join(syncfg.sailfishPath, 'lib'), '$LD_LIBRARY_PATH'])

outputDir = os.path.join(wd, prefix+'_quant')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

cmd = ' '.join(['sailfish quant -i', index, '-l "T=PE:O=><:S='+args.strand+'" -1', R1file, '-2', R2file, '-o', outputDir, '--threads 8'])
if not os.path.exists(os.path.join(wd, prefix+'_quant.sf')) and not os.path.exists(os.path.join(outputDir, 'quant.sf')):
	print 'sailfish start %s' % time.asctime()
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
	print 'sailfish end %s' % time.asctime()
	
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

# Load bias-corrected file to Synapse File.
print 'Loading %s to Synapse.' % os.path.join(wd, prefix+'_quant_bias_corrected.sf')
quantBCEntity = File(path=os.path.join(wd, prefix+'_quant_bias_corrected.sf'), name=prefix+'_quant_bias_corrected.sf', description='Quantified transcript isoforms.', parentId=args.out, synapseStore=True)	
quantBCEntity = syn.store(quantBCEntity, forceVersion=False, activity=act)
syn.setAnnotations(quantBCEntity, annotations=dict(fileType='quantitation',normalized='TPM',summaryLevel='transcript',biasCorrection='True'))
print 'new entity id %s' % quantBCEntity.id



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



# Move output files to head node for consolidation
evalName = 'QUANTsailfish'
hd = os.path.join(syncfg.headNFSPath, evalName)
if not os.path.exists(hd):
	os.mkdir(hd)
shutil.copy2(os.path.join(wd, prefix+'_quant.sf'), hd)
shutil.copy2(os.path.join(wd, prefix+'_quant_bias_corrected.sf'), hd)
shutil.copy2(os.path.join(wd, prefix+'_count_info.sf'), hd)


# clean up local files: BAM, FASTQ, sailfish files
os.remove(localBAMfilePath)
os.remove(R1file)
os.remove(R2file)
shutil.rmtree(wd)

# change status of BAM 
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)
