#! /usr/bin/env python
# Oct. 22, 2013
# KKD for Sage Bionetworks

## Hardcoded: hg19 genome, evalID
## Requires: varscan and samtools to be in the path
## To Do: return entity from add_workflow_step so that I don't have to fetch it again from synapse

import synapseclient, os, argparse, os.path, sys, boto
from synapseclient import *
sys.path.append('/home/kkdang/bin/synapseseq')
import seq_loading as sl

parser = argparse.ArgumentParser(description='Runs Synapse somatic SNV workflow using VarScan.')
parser.add_argument('--input-tumor', dest='tbam', required=True, help='Tumor Synapse BAM entity ID to process.')
parser.add_argument('--input-normal', dest='nbam', required=True, help='Normal Synapse BAM entity ID to process.')
parser.add_argument('--upload', dest='id', required=True, help='Synapse ID of folder to which results are uploaded .')
parser.add_argument('--local', dest='wd', required=False, help='Local directory for output files.', default=os.getcwd())
parser.add_argument('--bucket', dest='aws', required=False, help='If the data is in an external bucket, provide the name.', default=None)
args = parser.parse_args()

# varscan and samtools executable paths or load statements?

syn = synapseclient.Synapse()
syn.login()


def getFilePath(inBAMentity):
	if args.aws:
		s3 = boto.connect_s3()
		H3bucket = s3.get_bucket(args.aws) 

		keyName = inBAMentity.externalURL.lstrip('file:/')
		key = H3bucket.get_key(keyName)
		filePath = os.path.join(args.wd, os.path.basename(keyName))
		if not os.file.exists(filePath):
			key.get_contents_to_filename(filePath)	
	else:
		originalPath = BAMentity.externalURL
		(dir, sampleName) = os.path.split(originalPath)
		filePath = os.path.join('/work/DAT_116__ICGC-TCGA_seq-breakpoints_challenge/Data', sampleName) 

	print '%s' % filePath	
	return filePath


NBAMentity = syn.get(args.nbam, downloadFile = False)
normalFilePath = getFilePath(NBAMentity)
normalPrefix = os.path.basename(normalFilePath).rstrip('.bam')

TBAMentity = syn.get(args.tbam, downloadFile = False)
tumorFilePath = getFilePath(TBAMentity)
tumorPrefix = os.path.basename(tumorFilePath).rstrip('.bam')
	

### Run Varscan
# in: normalFilePath, tumorFilePath
# out: wigFile
wigFile = os.path.join(args.wd, prefix+'.bigwig')
bam_to_bigwig(bam=filePath, genome='hg19', output=wigFile)

normalPileup = ' '.join(['samtools mpileup -q 1 -f', reference, normalFilePath])
tumorPileup = ' '.join(['samtools mpileup -q 1 -f', reference, tumorFilePath])

cmd = ' '.join(['java -jar', varscan, 'somatic <', normalPileup, '<', tumorPileup, ' output'])

sys.exit()

### Load result to synapse
# params needed from input: pid for output, synid of this code

# specifics for this workflow step
provDict = dict()
provDict['folder'] = 'browserTracks'
provDict['softwareName'] = 'pybedtools'
provDict['description'] = 'Create coverage track from a BAM file in reads per million units.'
provDict['actName'] = 'generate coverage track'
provDict['fileDescription'] = 'bigWig-format file of read coverage in RPM units.'
provDict['executed'] = 'syn2286692' # syn2286692 is pybedtools
provDict['used'] = BAMentity
provDict['suffix'] = 'bigwig'
provDict['store'] = 'true'
annotations = dict(['fileType','bigwig'], ['normalized','RPM'])
provDict['annotations'] = annotations

browserCoverageEntityID = sl.add_workflow_step_to_synapse(wigFile, stepDict=provDict, parentid=args.id, syn=syn)


### Submit result to synapse BrowserCoverage Eval
browserCoverageEntity = syn.get(browserCoverageEntityID, downloadFile = False)
browserCoverageEval = syn.getEvaluation('2275609')
#tmp = syn.joinEvaluation(snpEval)
profile = syn.getUserProfile() 

submission = syn.submit(entity=browserCoverageEntity, evaluation = browserCoverageEval.id,  name = browserCoverageEntity.name, teamName = profile['displayName'])
print 'Submitted %s to %s' % (browserCoverageEntity.name, browserCoverageEval.name)
