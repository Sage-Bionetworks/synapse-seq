#! /usr/bin/env python
# Dec. 9, 2013
# KKD for Sage Bionetworks


import synapseclient, os, yaml, argparse
import seq_running as sr

parser = argparse.ArgumentParser(description='Checks for any new submissions to seq workflow evaluations and runs the corrsponding compute jobs.')
parser.add_argument('--local', dest='bin', required=False, help='Location for code files on this compute system.', default='~/bin/scripts')
args = parser.parse_args()

syn = synapseclient.Synapse()
syn.login()
profile = syn.getUserProfile() 
evalListFile = syn.get('syn2319132', downloadFile = True) # yaml file

stream = open(evalListFile.path, 'r')
for eval,code in yaml.load(stream).iteritems():
	print 'eval %s: code file: %s' % (eval, code)
	
	# Get code from local install
	codePath = os.path.join(args.bin, code)

	# Make location for logs 
	logsDir = os.path.join(os.getcwd(), 'logs')
	if not os.path.exists(logsDir):
		os.mkdir(logsDir)

	# Run count job for open submissions
	for submission in syn.getSubmissions(eval):	
#		print '%s' % submission.entityId
		submissionAnnotations = syn.getAnnotations(submission.entityId, version=submission.versionNumber)
		sr.runJobsForSubmissions(submission, codePath, logsDir, outputProjectID=submissionAnnotations['outputID'], commandLineParams=submissionAnnotations['commandLineParams'], externalBucket=submissionAnnotations['bucket'], syn=syn)


