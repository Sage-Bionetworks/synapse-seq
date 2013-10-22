#! /usr/bin/env python
# Oct. 18, 2013
# KKD for Sage Bionetworks

import synapseclient, os, argparse, sys, os.path, subprocess
from synapseclient import *

parser = argparse.ArgumentParser(description='Checks for new BAM submissions and runs them through synapse fusion evaluation.')
parser.add_argument('--BAMeval', dest='ID', required=True, help='Synapse ID of BAM evaluation that provides input files for the fusion evaluation.', default=None)
parser.add_argument('--bucket', dest='aws', required=False, help='If the data is in an external bucket, provide the name.', default=None)
args = parser.parse_args()

#### Global constants
syn = synapseclient.Synapse()
syn.login()
profile = syn.getUserProfile() 
bamEval = syn.getEvaluation(args.ID)


#### Process data	

# Get code from synapse for fusion eval
# evalCode = syn.get('syn999999', downloadFile = True, downloadFileLocation = /some/path)
evalCodePath = '/home/kristen/bin/synapseseq/scripts/eval_fusion_qsub.py'

# Get dict of existing submissions
for submission in syn.getSubmissions(bamEval.id):	

	status = syn.getSubmissionStatus(submission)
	if status.status == 'OPEN':

		if args.aws is not None:
			fusion_cmd = ' '.join(['qsub -N fusion -o /home/kristen/logs/$JOB_NAME.$JOB_ID -j y -S /usr/bin/python -V', evalCodePath, '--input', submission.entityId, '--params /home/kristen/ccle_fusion.config --config /home/kristen/fusionmap_generic.config --bucket', args.aws, '--local /mnt/kristen'])
		else:
			fusion_cmd = ' '.join(['qsub -N fusion -o /home/kristen/logs/$JOB_NAME.$JOB_ID -j y -S /usr/bin/python -V', evalCodePath, '--input', submission.entityId, '--params /home/kristen/ccle_fusion.config --config /home/kristen/fusionmap_generic.config --local /mnt/kristen'])

		print '%s' % fusion_cmd
		subprocess.call(fusion_cmd, shell = True)

# 		status.status = 'SCORED'
# 		status = syn.store(status)
				
