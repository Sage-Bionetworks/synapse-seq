#! /usr/bin/env python
# 29 Jul 2016
# Kristen K Dang for Sage Bionetworks


import os
import argparse
import subprocess
import time
import shutil
import yaml
import random
import sys
import synapseclient
from synapseclient import File, Table, Activity
from synapseseq import seq_running as sr
from toil.job import Job





# TODO Implement support for synapse queuing system, though don't require
def getSubmission():
	submission = syn.getSubmission(args.submission, downloadFile = False)
	localInputFilePath = os.path.join(wd, submission.name)
	if not os.path.exists(localInputFilePath):
		submission = syn.getSubmission(args.submission, downloadFile = True, downloadLocation = wd)
	prefix = os.path.basename(localInputFilePath).rstrip('.bam')
	cfPath = os.path.join(wd, '_'.join([prefix, evaluation.name, 'commands.txt']))
	commandsFile = open(os.path.join(wd, cfPath),'w')
	print >> commandsFile, '%s' % localInputFilePath		



### SamtoFastq
def samToFastq():
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



### Sailfish
def runSailfish():

	## Load output BAM file 
	print 'Loading %s to Synapse.' % outBAMfile
	outBAMEntity = File(path=outBAMfile, name=os.path.basename(outBAMfile), description='Aligned reads in BAM format.', parentId=config['star']['output'], synapseStore=True)
	outBAMEntity = syn.store(outBAMEntity, forceVersion=False)
	#syn.setAnnotations(outBAMEntity, annotations=config['star']['annotations'])
	print 'new entity id %s' % outBAMEntity.id


	## Star metrics
	table = syn.get(config['star']['metricsTable'])
	metrics = [prefix, args.submission]
	outLogFile = os.path.join(outputDir, prefix+'.Log.final.out')
	with open(outLogFile) as metricsFile:	
		for line in metricsFile:
			vals = line.split('|')
			if len(vals) > 1:
				metrics.append(vals[1].strip().rstrip('%'))
	metricsFile.close()
	syn.store(Table(table, [metrics]))


	## Set up provenance.
	commandsFile.close()
	print 'Loading %s to Synapse.' % cfPath
	cf = File(path=cfPath, description='Job commands.', parentId=config['workflow']['output'], synapseStore=True)
	cf = syn.store(cf)

	## Set provenance on BAM
	align = Activity(name='Alignment', description='Alignment of reads to human genome.', executed=[cf.id, 'syn2243150', 'syn4566096'])
	align.used({'reference':{'targetId':submission.entityId, 'targetVersionNumber':submission.versionNumber}})
	align.used(config['star']['ref'])
	syn.setProvenance(outBAMEntity, align)
 

	## clean up local files
	if 'keep-local' not in config['workflow']:
		os.remove(localInputFilePath)
		shutil.rmtree(wd)


def changeStatus():

	## change status of input submission
	status = syn.getSubmissionStatus(submission)
	status.status = 'SCORED' # Scored is functioning as "finished" for now.
	status = syn.store(status)



##################
# Workflow
##################


if __name__=="__main__":

	parser = argparse.ArgumentParser(description='Submits a workflow job to a grid queue manager such as SGE.')
	parser.add_argument('submission', help='Submission id of input file.')
	parser.add_argument('--eval', required=True, help='Synapse ID of the eval to run.')
	parser.add_argument('--conf', required=True, help='Path to local config file.')
	parser.add_argument('--debug', required=False, help='Print all commands to output file even if not executed.', action='store_true')
	args = parser.parse_args()


	## Load workflow params
	config = yaml.load(file(args.conf))


	## Log in to Synapse
	syn = synapseclient.login()
	evaluation = syn.getEvaluation(args.eval)



	## Wrap jobs
	
	
	
	## specify order
	
	Job.Runner.startToil(j1, args)