#! /usr/bin/env python
# Oct. 17, 2013
# KKD for Sage Bionetworks

import synapseclient, os, argparse, subprocess, os.path, sys
from synapseclient import *
sys.path.append('/Users/kristen/Computing/synapseSeq/')
import seq_loading as sl

parser = argparse.ArgumentParser(description='Runs Synapse fusion reads prediction workflow using omicsoft FusionMap.')
parser.add_argument('--input', dest='bam', required=True, help='Synapse BAM entity ID to process.')
parser.add_argument('--params', dest='path', required=True, help='Path to config file for synapse provenance and annotations specific to this project.')
parser.add_argument('--local', dest='wd', required=False, help='Local directory for output files.', default=os.getcwd())
parser.add_argument('--config', dest='cfg', required=False, help='Path to default FusionMap config file.', default=os.path.join(os.getcwd(), 'generic.config'))
args = parser.parse_args()

fusionExe = 'mono ~/Software/FusionMap_2013-07-30/bin/FusionMap.exe'
fusionDir = '~/Software/FusionMap_2013-07-30/'

syn = synapseclient.Synapse()
syn.login()


BAMentity = syn.get(args.bam, downloadFile = False)
if 'path' in BAMentity:
	filePath = BAMentity.path
else:
	filePath = BAMentity.name
	
prefix = os.path.basename(filePath).rstrip('.bam')

### Get unmapped reads
tmpfile = os.path.join(args.wd, prefix+'_tmp1.bam')
cmd = ' '.join(['samtools view -u -f 4 -F264', filePath, '>', tmpfile])
print '%s' % cmd
tmpfile = os.path.join(args.wd, prefix+'_tmp2.bam')
cmd = ' '.join(['samtools view -u -f 8 -F260', filePath, '>', tmpfile])
print '%s' % cmd
tmpfile = os.path.join(args.wd, prefix+'_tmp3.bam')
cmd = ' '.join(['samtools view -u -f 12 -F256', filePath, '>', tmpfile])
print '%s' % cmd

cmd = ''.join(['samtools merge -u - ', args.wd, '/', prefix+'tmp[123].bam | samtools sort -n - unmapped'])
print '%s' % cmd


### Run FusionMap on unmapped reads

# Make fusionmap config file
configFile = open(os.path.join(args.wd, '.'.join([prefix, 'config'])), 'w')
configFile.write('\n'.join(['<Files>', filePath, '\n']))
generic = open(args.cfg, 'r')
for line in generic:
	configFile.write(line)
configFile.write('\n'.join(['', '<Output>', 'OutputPath='+args.wd, 'OutputName='+prefix]))
configFile.close()
generic.close()

# Run FusionMap
cmd = ' '.join([fusionExe, '--semap', fusionDir, 'GRCh37 Gencode14', os.path.join(args.wd, '.'.join([prefix, 'config']))])
print '%s' % cmd


### Load result to synapse
loadFilePath = os.path.join(args.wd, prefix+'.FusionReport.txt')
print '%s' % loadFilePath

synDict = sl.readConfig(args.path)
fullUsed = ','.join([synDict['1']['used'], BAMentity.id])
synDict['1']['used'] = fullUsed
fusionEntityID = sl.add_workflow_step_to_synapse(loadFilePath, stepDict=synDict['1'], parentid=synDict['projectID'], syn=syn)


### Submit result to synapse Fusion Eval
fusionEntity = syn.get(fusionEntityID, downloadFile = False)
fusionEval = syn.getEvaluation('2275609')
#tmp = syn.joinEvaluation(snpEval)
profile = syn.getUserProfile() 

submission = syn.submit(entity=fusionEntity, evaluation = fusionEval.id,  name = fusionEntity.name, teamName = profile['displayName'])
print 'Submitted %s to %s' % (fusionEntity.name, fusionEval.name)
