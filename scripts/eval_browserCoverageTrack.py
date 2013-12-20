#! /usr/bin/env python
# Oct. 22, 2013
# KKD for Sage Bionetworks

## Hardcoded: hg19 genome, evalID
## Requires: bedtools and samtools to be in the path
## To Do: return entity from add_workflow_step so that I don't have to fetch it again from synapse

import synapseclient, os, argparse, os.path, sys, boto
from pybedtools.contrib.bigwig import bam_to_bigwig
from synapseclient import *
sys.path.append('/home/kkdang/bin/synapseseq')
import seq_loading as sl

parser = argparse.ArgumentParser(description='Runs Synapse genome browser coverage track creation workflow using pybedtools.')
parser.add_argument('--input', dest='bam', required=True, help='Synapse BAM entity ID to process.')
parser.add_argument('--upload', dest='id', required=True, help='Synapse ID of folder to which results are uploaded .')
parser.add_argument('--local', dest='wd', required=False, help='Local directory for output files.', default=os.getcwd())
parser.add_argument('--bucket', dest='aws', required=False, help='If the data is in an external bucket, provide the name.', default=None)
args = parser.parse_args()

# bedtools and samtools executable paths or load statements?

syn = synapseclient.Synapse()
syn.login()

BAMentity = syn.get(args.bam, downloadFile = False)

if args.aws:
	s3 = boto.connect_s3()
	H3bucket = s3.get_bucket(args.aws) 

	keyName = BAMentity.externalURL.lstrip('file:/')
	key = H3bucket.get_key(keyName)
	filePath = os.path.join(args.wd, os.path.basename(keyName))
	if not os.file.exists(filePath):
		key.get_contents_to_filename(filePath)	
else:
	originalPath = BAMentity.externalURL
	(dir, fileName) = os.path.split(originalPath)
	(sampleDir, subdir) = os.path.split(dir)
	(rootDir, sampleName) = os.path.split(sampleDir)
	filePath = os.path.join('/work/DAT_112__Craniosynostosis_seq_subtype/Staging/In/sage', sampleName, subdir, fileName) 

print '%s' % filePath	
prefix = os.path.basename(filePath).rstrip('.bam')

### Make bigWig file
# in: filePath
# out: wigFile
wigFile = os.path.join(args.wd, prefix+'.bigwig')
bam_to_bigwig(bam=filePath, genome='hg19', output=wigFile)


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
