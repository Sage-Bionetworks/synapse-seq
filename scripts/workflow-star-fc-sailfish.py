#! /usr/bin/env python
# June 2, 2015
# Kristen K Dang for Sage Bionetworks
# Requires: UrQt, star, picard, featurecounts
# Uses max 31GB RAM


import os, argparse, subprocess, sys, time, shutil, yaml, random
import synapseclient, synapseseq
from synapseclient import File, Table, Activity
from synapseseq import seq_running as sr
from picard_parser import PicardParser



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
localInputFilePath = os.path.join(wd, submission.name)
if not os.path.exists(localInputFilePath):
	submission = syn.getSubmission(args.submission, downloadFile = True, downloadLocation = wd)
prefix = os.path.basename(localInputFilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evaluation.name, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print >> commandsFile, '%s' % localInputFilePath		



## Run workflow

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



### Quality trim FASTQ files
random.seed()
tmpPre = str(random.randint(0, 9999999))

if config['workflow']['paired'] is False:
	R1TrimFile = os.path.join(wd, prefix + '_SE_trimmed.fq.gz')
	tempR1TrimFile = os.path.join(wd, tmpPre + '_SE.fq.gz')
	trimOutFile = os.path.join(wd, prefix + '_trim.stdout')
	R2TrimFile = ''

	cmd = ' '.join(['UrQt --in', R1file, '--out', tempR1TrimFile, '--gz --m', str(config['workflow']['threads']), '--min_read_size 30 --phred 33 --t 20 --buffer 100000 >', trimOutFile])

	if not os.path.exists(R1TrimFile):
		print >> commandsFile, '%s' % cmd
		subprocess.call(cmd, shell = True)
		os.rename(tempR1TrimFile, R1TrimFile)
	elif args.debug is True:
		print >> commandsFile, '%s' % cmd	

else:
	R1TrimFile = os.path.join(wd, prefix + '_R1_trimmed.fq.gz')
	R2TrimFile = os.path.join(wd, prefix + '_R2_trimmed.fq.gz')
	tempR1TrimFile = os.path.join(wd, tmpPre + '_R1.fq.gz')
	tempR2TrimFile = os.path.join(wd, tmpPre + '_R2.fq.gz')
	trimOutFile = os.path.join(wd, prefix + '_trim.stdout')

	cmd = ' '.join(['UrQt --in', R1file, '--inpair', R2file, '--out', tempR1TrimFile, '--outpair', tempR2TrimFile, '--gz --m', str(config['workflow']['threads']), '--min_read_size 30 --phred 33 --t 20 --buffer 100000 >', trimOutFile])

	if not os.path.exists(R1TrimFile) or not os.path.exists(R2TrimFile):
		print >> commandsFile, '%s' % cmd
		subprocess.call(cmd, shell = True)
		os.rename(tempR1TrimFile, R1TrimFile)
		os.rename(tempR2TrimFile, R2TrimFile)
	elif args.debug is True:
		print >> commandsFile, '%s' % cmd	
	

## Trim metrics
table = syn.get(config['urqt']['metricsTable'])
metrics = [prefix, args.submission]
with open(trimOutFile) as metricsFile:	
	for line in metricsFile:
		vals = line.split(':')
		if len(vals) == 2:
			data_vals = vals[1].strip().split()
			if len(data_vals) == 1:
				metrics.append(vals[1].strip())
			if len(data_vals) == 2:
				metrics.append(data_vals[1].strip('()%'))
		if config['workflow']['paired'] is False:
			metrics.extend(list('NA', 'NA', 'NA', 'NA', 'NA'))
metricsFile.close()
syn.store(Table(table, [metrics]))



### Run STAR
outputDir = os.path.join(wd, prefix+'_align')
if not os.path.exists(outputDir):
	os.mkdir(outputDir) 

index = os.path.join(config['system']['headNFSPath'], config['star']['index'])

cmd = ' '.join(['STAR --runMode alignReads --runThreadN', str(config['workflow']['threads']), '--genomeDir', index, '--readFilesIn', R1TrimFile, R2TrimFile, '--readFilesCommand zcat --outFileNamePrefix', os.path.join(outputDir,prefix+'.'), '--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within'])
outBAMfile = os.path.join(outputDir, prefix+'.Aligned.sortedByCoord.out.bam')
if not os.path.exists(outBAMfile):
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
elif args.debug is True:
	print >> commandsFile, '%s' % cmd


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




### Run Picard metrics
#For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand. Required. Possible values: {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}
if config['workflow']['stranded'] is 'First':
	strand = 'FIRST_READ_TRANSCRIPTION_STRAND'
elif config['workflow']['stranded'] is 'Second':
	strand = 'SECOND_READ_TRANSCRIPTION_STRAND'
else:
	strand = 'NONE'

cmd = ' '.join(['java -Xmx2G -jar', os.path.join(config['system']['picard'], 'CollectRnaSeqMetrics.jar'), 'REF_FLAT=', os.path.join(config['system']['headNFSPath'], config['picardRNA']['refflat']), ' RIBOSOMAL_INTERVALS=', os.path.join(config['system']['headNFSPath'], config['picardRNA']['rRNAintervals']), 'STRAND_SPECIFICITY='+strand, 'CHART_OUTPUT=', os.path.join(wd, os.path.basename(outBAMfile)+'.pdf'), 'INPUT=', outBAMfile, 'OUTPUT=', os.path.join(wd, os.path.basename(outBAMfile)+'.picardRNA'), 'ASSUME_SORTED=true'])

if not os.path.exists(os.path.join(wd, os.path.basename(outBAMfile)+'.picardRNA')):
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
elif args.debug is True:
	print >> commandsFile, '%s' % cmd

## PicardRNA
pp  = PicardParser(filename=os.path.join(wd,os.path.basename(outBAMfile)+'.picardRNA'))
table = syn.get(config['picardRNA']['metricsTable'])
metrics = [prefix, args.submission]
for item in pp.metrics.itervalues():
	metrics.append(item)
syn.store(Table(table, [metrics]))




### Run featurecounts
outCountsFile = os.path.join(wd, prefix+'_gene_counts.txt')
gtf = os.path.join(config['system']['headNFSPath'], config['featurecounts']['gtf'])
os.chdir(wd) # Attempt to redirect temp output files to this location, since they seem to go to the wd with no option to specify elsewhere.

# strand specificity: 
# Indicate if strand-specific read counting should be performed. It has three possible values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default. For paired-end reads, strand of the first read is taken as the strand of the whole fragment and FLAG field of the current read is used to tell if it is the first read in the fragment
if config['workflow']['stranded'] is 'First':
	strand = '1'
elif config['workflow']['stranded'] is 'Second':
	strand = '2'
else:
	strand = '0'

# paired end
if config['workflow']['paired'] is False:
	fc = 'featureCounts'
else:
	fc = 'featureCounts -p'

# exon vs gene		
if (config['featurecounts']['count-unit'] is 'exon'): # If exon counting specified
	cmd = ' '.join([fc, '-t exon -g gene_id -a', gtf, '-o', outCountsFile, '-s', strand, '-T', config['workflow']['threads'], '-f', outBAMfile])
else:
	cmd = ' '.join([fc, '-t exon -g gene_id -a', gtf, '-o', outCountsFile, '-s', strand, '-T', str(config['workflow']['threads']), outBAMfile])

if not os.path.exists(outCountsFile):
	print >> commandsFile, '%s' % cmd
	subprocess.call(cmd, shell = True)
elif args.debug is True:
	print >> commandsFile, '%s' % cmd

## featurecounts
table = syn.get(config['featurecounts']['metricsTable'])
metrics = [prefix, args.submission]
outLogFile = os.path.join(wd, prefix+'_gene_counts.txt.summary')
with open(outLogFile) as metricsFile:	
	for line in metricsFile:
		if line.startswith('Status'): continue
		metrics.append(line.split()[1])
metricsFile.close()
syn.store(Table(table, [metrics]))



### Run sailfish
# ********

	
commandsFile.close()


# Load results to synapse

## Set up provenance.
print 'Loading %s to Synapse.' % cfPath
cf = File(path=cfPath, description='Job commands.', parentId=config['workflow']['output'], synapseStore=True)
cf = syn.store(cf)

## Set provenance on BAM
align = Activity(name='Alignment', description='Alignment of reads to human genome.', executed=[cf.id, 'syn2243150', 'syn4566096'])
align.used({'reference':{'targetId':submission.entityId, 'targetVersionNumber':submission.versionNumber}})
align.used(config['star']['ref'])
syn.setProvenance(outBAMEntity, align)

## Load output counts file 
print 'Loading %s to Synapse.' % outCountsFile
outCountsEntity = File(path=outCountsFile, name=os.path.basename(outCountsFile), description='Gene or exon counts based on alignment to human genome.', parentId=config['featurecounts']['output'], synapseStore=True)	
outCountsEntity = syn.store(outCountsEntity, forceVersion=False, activityName='Quantitation', activityDescript='Counting reads that align to gene models.', used=[config['featurecounts']['ref']], executed=[cf.id, 'syn2807330'])
outCountsEntity = syn.store(outCountsEntity, forceVersion=False, activityName='Quantitation', activityDescript='Counting reads that align to gene models.', used=[outBAMEntity, config['featurecounts']['ref']], executed=[cf.id, 'syn2807330'])
#syn.setAnnotations(outCountsEntity, annotations=config['featurecounts']['annotations'])
print 'new entity id %s' % outCountsEntity.id

 

## clean up local files
if 'keep-local' not in config['workflow']:
	os.remove(localInputFilePath)
	shutil.rmtree(wd)

## change status of input submission
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)

