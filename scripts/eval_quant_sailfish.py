#! /usr/bin/env python
# Sept. 26, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: sailfish, picard
# Uses ~ 5GB RAM


import os, argparse, subprocess, sys, time, shutil, yaml, random
import synapseclient, synapseseq
from synapseclient import Activity, File, Table
from synapseseq import seq_running as sr


parser = argparse.ArgumentParser(description='Quantifies transcripts in input BAM using sailfish.')
parser.add_argument('submission', dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--index', dest='idx', required=True, help='Synapse ID of transcript index required by sailfish.', default=None)
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
localInputFilePath = os.path.join(wd, submission.name)
if not os.path.exists(localInputFilePath):
	submission = syn.getSubmission(args.submission, downloadFile = True, downloadLocation = wd)
prefix = os.path.basename(localInputFilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evaluation.name, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localInputFilePath		



## Run workflow

## Get index
index = os.path.join(config['system']['headNFSPath'], config['sailfish']['index'])


### SamtoFastq
if config['workflow']['paired'] is False:
	R1file = os.path.join(wd, prefix + '_SE.fastq')
	cmd = ' '.join(['java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SortSam.jar'), 'INPUT=', localInputFilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'TMP_DIR=', os.path.join(wd, 'tmp'), 'VALIDATION_STRINGENCY=SILENT'])

	if not os.path.exists(R1file):
		print >> commandsFile, '%s' % cmd
		subprocess.call(cmd, shell = True)
	elif args.debug is True:
		print >> commandsFile, '%s' % cmd
else:
	R1file = os.path.join(wd, prefix + '_R1.fastq')
	R2file = os.path.join(wd, prefix + '_R2.fastq')
	cmd = ' '.join(['java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SortSam.jar'), 'INPUT=', localInputFilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file, 'TMP_DIR=', os.path.join(wd, 'tmp'), 'VALIDATION_STRINGENCY=SILENT'])

	if not os.path.exists(R1file) or not os.path.exists(R2file):
		print >> commandsFile, '%s' % cmd
		subprocess.call(cmd, shell = True)
	elif args.debug is True:
		print >> commandsFile, '%s' % cmd


## Run SAILFISH step
os.environ["PATH"] = ':'.join([os.path.join(config['system']['sailfishPath'], 'bin'), '$PATH'])
os.environ["LD_LIBRARY_PATH"] = ':'.join([os.path.join(config['system']['sailfishPath'], 'lib'), '$LD_LIBRARY_PATH'])

outputDir = os.path.join(wd, prefix+'_quant')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

cmd = ' '.join(['sailfish quant -i', index, config['sailfish']['params'], '-1', R1file, '-2', R2file, '-o', outputDir])
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
cf = File(path=cfPath, description='Job commands.', parentId=config['sailfish']['output'], synapseStore=True)
cf = syn.store(cf, activityName='quant_evaluation', executed=['https://github.com/Sage-Bionetworks/synapse-seq/blob/master/scripts/eval_quant_sailfish.py'])
act = Activity(name='transcript quantitation', description='Alignment-free transcript quantitation using Sailfish.', executed=['syn2325155', cf.id])
act.used(target=submission.entityId, targetVersion=submission.versionNumber)
act.used(config['sailfish']['ref'])

# Load raw quant file
print 'Loading %s to Synapse.' % os.path.join(wd, prefix+'_quant.sf')
quantEntity = File(path=os.path.join(wd, prefix+'_quant.sf'), name=prefix+'_quant.sf', description='Quantified transcript isoforms.', parentId=config['sailfish']['output'], synapseStore=True)	
quantEntity = syn.store(quantEntity, forceVersion=False, activity=act)
syn.setAnnotations(config['sailfish']['annotations'])
print 'new entity id %s' % quantEntity.id

# Load bias-corrected file to Synapse File.
print 'Loading %s to Synapse.' % os.path.join(wd, prefix+'_quant_bias_corrected.sf')
quantBCEntity = File(path=os.path.join(wd, prefix+'_quant_bias_corrected.sf'), name=prefix+'_quant_bias_corrected.sf', description='Quantified transcript isoforms.', parentId=config['sailfish']['output'], synapseStore=True)	
quantBCEntity = syn.store(quantBCEntity, forceVersion=False, activity=act)
syn.setAnnotations(config['sailfish']['annotations'])
print 'new entity id %s' % quantBCEntity.id



## Load metrics to Synapse table
table = syn.get(config['sailfish']['metricsTable'])
metrics = [prefix, args.submission]
outLogFile = os.path.join(wd, prefix+'_count_info.sf')
with open(outLogFile) as metricsFile:	
	for line in metricsFile:
		if line.startswith('Status'): continue
		metrics.append(line.split()[1])
metricsFile.close()
syn.store(Table(table, [metrics]))



# Move output files to head node for consolidation
evalName = 'QUANTsailfish'
hd = os.path.join(config['system']['headNFSPath'], evalName)
if not os.path.exists(hd):
	os.mkdir(hd)
shutil.copy2(os.path.join(wd, prefix+'_quant.sf'), hd)
shutil.copy2(os.path.join(wd, prefix+'_quant_bias_corrected.sf'), hd)
shutil.copy2(os.path.join(wd, prefix+'_count_info.sf'), hd)


## clean up local files
if 'keep-local' not in config['workflow']:
	os.remove(localInputFilePath)
	os.remove(R1file)
	os.remove(R2file)
	shutil.rmtree(wd)


## change status of input submission
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)
