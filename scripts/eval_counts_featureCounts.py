#! /usr/bin/env python
# Nov. 4, 2014
# Kristen K Dang for Sage Bionetworks
# Requires: featurecounts
# Uses ~ 60MB RAM, 4xinput BAM GB in disk


import os, argparse, subprocess, sys, time, shutil
import synapseclient, synapseseq
from synapseclient import Activity, File
from synapseseq import seq_running as sr
from synapseseq import synseqConfig as syncfg


parser = argparse.ArgumentParser(description='Counts reads aligned within boundaries of  genes or exons in input BAM using featurecounts.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--gtf', dest='gtf', required=True, help='Synapse ID of GTF defining gene models to use for counting.', default=None)
parser.add_argument('--output', dest='out', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--stranded', dest='strand', required=False, help='Strandedness of the reads (see featurecounts docs for explanation of options).', default='0')
parser.add_argument('--threads', dest='thread', required=False, help='Number of threads to use in featurecounts.', default='1')
parser.add_argument('--exon', dest='exon', required=False, help='Count exons rather than genes.', default=False)
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


subprocess.call(' '.join(['sudo chown ubuntu:ubuntu', syncfg.localWDPath]), shell = True)
sr.getSynConfigFromHead(localPath='/home/'+syncfg.workerUserName+'/',headPath=syncfg.headNFSPath)

## Set up output directory, log in to synapse	
evalName = 'COUNTfeaturecounts'
wd = os.path.join(syncfg.localWDPath, evalName)
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
gtf = os.path.join(syncfg.headNFSPath, gtfName)


## Get submission
submission = syn.getSubmission(args.bam, downloadFile = False)
localBAMfilePath = sr.getBAMtoComputeNode(wd=wd,syn=syn,submission=submission,bucket=args.bucket,extKey=args.key)
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')
cfPath = os.path.join(wd, '_'.join([prefix, evalName, 'commands.txt']))
commandsFile = open(os.path.join(wd, cfPath),'w')
print '%s' % localBAMfilePath		



## Run featurecounts step
outputFile = os.path.join(wd, prefix+'_gene_counts.txt')

os.chdir(wd) # Attempt to redirect temp output files to this location, since they seem to go to the wd with no option to specify elsewhere.
if (args.exon): # If exon counting specified
	cmd = ' '.join(['featureCounts -p -t exon -g gene_id -a', gtf, '-o', outputFile, '-s', args.strand, '-T', args.thread, '-f', localBAMfilePath])
else:
	cmd = ' '.join(['featureCounts -p -t exon -g gene_id -a', gtf, '-o', outputFile, '-s', args.strand, '-T', args.thread, localBAMfilePath])
print 'featurecounts start %s' % time.asctime()
print >> commandsFile, '%s' % cmd
subprocess.call(cmd, shell = True)
print 'featurecounts end %s' % time.asctime()
	


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



## Use this code after related JIRA is resolved
## Load metrics to Synapse table
# table = syn.get(syncfg.featurecountsMetricsTable)
# metrics = list()
# metrics.append(prefix)
# with open(os.path.join(wd, prefix+'_count_info.sf')) as metricsFile:	
# 	for line in metricsFile:
#		if not line.startswith('Status'):
#	 		metrics.append(line.split()[1])
# metricsFile.closed
# syn.store(Table(table, [metrics]))


# Move output files to head node for consolidation ... change this to dumping read counts into a table when Tables is ready
evalName = 'COUNTfeaturecounts'
hd = os.path.join(syncfg.headNFSPath, evalName)
if not os.path.exists(hd):
	os.mkdir(hd)
shutil.copy2(outputFile, hd)
shutil.copy2(outputFile+'.summary', hd)


# clean up BAM files
os.remove(localBAMfilePath)
os.remove(cfPath)
shutil.rmtree(wd)

# change status of BAM 
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "finished" for now.
status = syn.store(status)
