#! /usr/bin/env python
# Dec. 9, 2013
# KKD for Sage Bionetworks

# TODO Hardcoded: path to python in qsub statement.
# TODO Generalize qsub commands so not dependent on SGE, but e.g. PBS could also be used.
## User should manually move reference onto machine and verify integrity. The reference file should be in a location that is currently hardcoded but can be altered with a config file. The reference file should have the same name as the filename of the synapse entity or the basename of the external URL.

import synapseclient, os, yaml, argparse, subprocess

parser = argparse.ArgumentParser(description='Checks for any new submissions to seq workflow evaluations and runs the corrsponding compute jobs.')
parser.add_argument('--evals', dest='synid', required=False, help='Synapse ID of yaml-format file containing evaluation ids and corresponding script names.')
parser.add_argument('--synapseseq', dest='bin', required=False, help='Path to installation location of scripts files from synapseseq package.', default='/usr/local/somewhere')
args = parser.parse_args()

syn = synapseclient.login()
profile = syn.getUserProfile() 
evalListFile = syn.get(args.synid, downloadFile = True) # yaml file


def parseJobParameters(synapseID):
	yamlFileEntity = syn.get(synapseID)
	stream = open(yamlFileEntity.path, 'r')
	
	jobParams = ''
	for params in yaml.load(stream).iteritems():
		for item in params:
			jobParams += ' '+str(item)
	jobParams += ' --params '+synapseID
	print 'Read these params: %s' % jobParams
	return(jobParams)

paramFiles = dict()

stream = open(evalListFile.path, 'r')
for eval,code in yaml.load(stream).iteritems():
	print 'eval %s: code file: %s' % (eval, code)
	
	codePath = os.path.join(args.bin, code)

	logsDir = os.path.join(os.getcwd(), 'logs')
	if not os.path.exists(logsDir):
		os.mkdir(logsDir)

	# Run jobs for open submissions
	for submission in syn.getSubmissions(eval, status='RECEIVED'):	
#		print '%s' % submission.entityId
			additionalParams = ''
			entityAnnotations = syn.getAnnotations(submission.entityId)
			if 'bucket' in entityAnnotations:
				additionalParams += '--bucket ' + entityAnnotations['bucket'][0]
				additionalParams += ' --keyname ' + entityAnnotations['key'][0]
			if 'eval_'+str(eval) in entityAnnotations: # Are there command line params in the entity's annotations?
				evalName = 'eval_'+str(eval)
				if entityAnnotations[evalName][0] not in paramFiles: # check whether this file has been parsed
					additionalParams = parseJobParameters(entityAnnotations[evalName][0])
					paramFiles[entityAnnotations[evalName][0]] = additionalParams
				else:
					additionalParams = paramFiles[entityAnnotations[evalName][0]]

			cmd = ' '.join(['qsub -N synapseseq -o', os.path.join(logsDir, '$JOB_NAME.$JOB_ID'), '-j y -S /usr/bin/python -V', codePath, '--submission', submission.id, additionalParams])

			print '%s' % cmd
 			subprocess.call(cmd, shell = True)
# 			status.status = 'OPEN' 
# 			status = syn.store(status)



