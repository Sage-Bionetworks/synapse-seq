#! /usr/bin/env python
# Nov. 26, 2013
# KKD for Sage Bionetworks
# Requires: htseq installed globally

import synapseclient, os, argparse, subprocess, sys, boto, synapseseq
import seq_loading as sl

parser = argparse.ArgumentParser(description='Runs Synapse read counting workflow using HTSeq.')
parser.add_argument('--input', dest='bam', required=True, help='Synapse BAM entity ID to process.')
parser.add_argument('--output', dest='sid', required=True, help='Synapse ID of project or folder to contain the output data.')
parser.add_argument('--GTF', dest='gtf', required=False, help='Synapse ID of GTF file to use for gene models.', default='syn2315951') # syn2315951 is ENSGv72
parser.add_argument('--local', dest='wd', required=False, help='Local directory for output files.', default=os.getcwd())
parser.add_argument('--count-mode', dest='count', required=False, help='Counting mode for HTSeq.', default='intersection-nonempty')
parser.add_argument('--stranded', dest='strand', required=False, help='Stranded counting in HTSeq?', default='no')
args = parser.parse_args()


syn = synapseclient.Synapse()
syn.login()

BAMentity = syn.getSubmission(args.bam, downloadFile = False)
BAMannotations = syn.getAnnotations(BAMentity.entityId)

GTFentity = syn.get(args.gtf, downloadFile = False)
print 'Using gene models found in %s' % GTFentity.name

if ('externalURL' in GTFentity) and (GTFentity.externalURL.startswith('ftp')):
	cmd = ' '.join(['wget', GTFentity.externalURL])
	print '%s' % cmd
	subprocess.call(cmd, shell = True)
	GTFpath = os.path.join(args.wd, os.path.basename(GTFentity.externalURL))
	if GTFentity.externalURL.endswith('gz'):
		subprocess.call(' '.join(['gunzip', os.path.basename(GTFentity.externalURL)]), shell = True)
		GTFpath = os.path.join(args.wd, os.path.basename(GTFentity.externalURL).rstrip('.gz'))
else:
	GTFentity = syn.get(args.gtf, downloadFile = True, downloadLocation = args.wd)
	GTFpath = os.path.join(args.wd, GTFentity.name)	

if 'bucket' in BAMannotations:
	s3 = boto.connect_s3()
	H3bucket = s3.get_bucket(BAMannotations['bucket'][0]) 

	keyName = BAMannotations['key'][0]
	bucketItem = H3bucket.get_key(keyName)
	filePath = os.path.join(args.wd, os.path.basename(keyName))
	if not os.path.exists(filePath):
		print 'Getting data in bucket %s' % BAMannotations['bucket'][0]
		bucketItem.get_contents_to_filename(filePath)	
else:
	filePath = os.path.join(args.wd, BAMentity.name)
	if not os.path.exists(filePath):
		print 'Will download %s from synapse.' % BAMentity.name
		BAMentity = syn.get(args.bam, downloadFile = True, downloadLocation = args.wd)

print '%s' % filePath	
prefix = os.path.basename(filePath).rstrip('.bam')

provDict = dict()
provDict['folder'] = 'counts'
provDict['softwareName'] = 'HTSeq'
provDict['description'] = 'Counting reads aligned to gene models.'
provDict['actName'] = 'Counting'
provDict['fileDescription'] = 'Raw read counts.'
provDict['executed'] = 'syn2243147,1' # syn2243147 is htseq
provDict['suffix'] = 'htseq'
provDict['store'] = 'true'
provDict['used'] = ','.join([BAMentity.id, GTFentity.id])
provDict['annotations'] = dict(fileType='count',normalized='no',summaryLevel='gene')


### Run samsort step if necessary
if ('sorted' not in BAMannotations) or (BAMannotations['sorted'][0] != 'queryname'):
	filePath = ''.join([filePath.rstrip('.bam'), '_namesort.bam']) 
	if not os.path.exists(filePath):
		cmd = ' '.join(['samtools sort -n', filePath, prefix+'_namesort'])
		print '%s' % cmd
		subprocess.call(cmd, shell = True)
	
### Run htseq step
outputFile = filePath.rstrip('.bam') + '.htseq'				
if not os.path.exists(outputFile):
	cmd = ' '.join(['samtools view ', filePath, '| python -m HTSeq.scripts.count -m', args.count, '-s', args.strand, '-', GTFpath, '>', outputFile])
	print '%s' % cmd
	subprocess.call(cmd, shell = True)


### Load result to synapse
print '%s' % outputFile
countEntityID = sl.add_workflow_step_to_synapse(os.path.join(args.wd, outputFile), stepDict=provDict, parentid=args.sid, syn=syn)

status = syn.getSubmissionStatus(BAMentity)
status.status = 'CLOSED' # Closed is functioning as "completed" for now.
status = syn.store(status)

### Submit result to synapse counted Eval?
# countEntity = syn.get(countEntityID, downloadFile = False)
# countEval = syn.getEvaluation('2313942') # 2313942 is counts_HTseq
# profile = syn.getUserProfile() 
# submission = syn.submit(entity=countEntity, evaluation = countEval.id,  name = countEntity.name, teamName = profile['displayName'])
# print 'Submitted %s to %s' % (countEntity.name, countEval.name)
