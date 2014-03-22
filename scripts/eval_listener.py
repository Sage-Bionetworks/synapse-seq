#! /usr/bin/env python
# Dec. 9, 2013
# KKD for Sage Bionetworks
# Requires: synapseseq installed such that contents of scripts folder are callable on command line
# To do: record version of synapseseq in provenance

import synapseclient, os, yaml, argparse, synapseseq
import seq_running as sr

parser = argparse.ArgumentParser(description='Checks for any new submissions to seq workflow evaluations and runs the corrsponding compute jobs.')
parser.add_argument('--evals', dest='synid', required=False, help='Synapse ID of yaml-format file containing evaluation ids and corresponding script names.', default='syn2319132')
args = parser.parse_args()

syn = synapseclient.Synapse()
syn.login()
profile = syn.getUserProfile() 
evalListFile = syn.get(args.synid, downloadFile = True) # yaml file

stream = open(evalListFile.path, 'r')
for eval,code in yaml.load(stream).iteritems():
	print 'Checking submissions for eval %s: %s' % (eval, code)
	
	# Make location for logs 
	logsDir = os.path.join(os.getcwd(), 'logs')
	if not os.path.exists(logsDir):
		os.mkdir(logsDir)

	# Run count job for open submissions
	for submission in syn.getSubmissions(eval):	
#		print '%s' % submission.entityId
		submissionAnnotations = syn.getAnnotations(submission.entityId, version=submission.versionNumber)
		
		if 'bucket' in submissionAnnotations:
			sr.runJobsForSubmissions(submission, code, logsDir, outputProjectID=submissionAnnotations['outputID'], commandLineParams=submissionAnnotations['workflowParams'], externalBucket=submissionAnnotations['bucket'], syn=syn)

		else:
			sr.runJobsForSubmissions(submission, code, logsDir, outputProjectID=submissionAnnotations['outputID'], commandLineParams=submissionAnnotations['workflowParams'], syn=syn)
