"""
Functions for running sequencing workflows using Synapse.

"""

# To do items: 
#	unit tests

import synapseclient, os, argparse, sys, os.path, subprocess
import seq_loading as sl


def getBAMs(projectOrFolderID, syn):
	'''Returns list of BAM entities in the given project or folder id.'''
	
	# Check if projectOrFolderID is BAM folder
	container = syn.get(projectOrFolderID, downloadFile = False)
	if container.name == 'BAM':
		BAMfolderID = container.id	
	# If not, get ID of BAM folder for this container
	else:
		allFoldersDict = sl.getExistingFolders(syn,container.id)
		BAMfolderID = allFoldersDict['BAM']
	
	BAMEntityList = list()
	print 'BAM folder ID: %s' % BAMfolderID
	sql = ''.join(['SELECT id, name FROM entity WHERE parentId=="', BAMfolderID, '"'])
	print '%s' % sql
	results = syn.chunkedQuery(sql)
	for result in results:
		if result['entity.name'].endswith('.bam'):
			BAMEntityList.append(syn.get(result['entity.id'], downloadFile = False))
	
	return(BAMEntityList)
	


def runJobsForSubmissions(submission, evalCode, logsDir, outputProjectID, commandLineParams, syn, externalBucket=None):
	'''	Checks for new submissions and runs corresponding evaluation jobs. Requires qsub.'''
	# Hardcoded: path to python in qsub statement.


	profile = syn.getUserProfile() 
	status = syn.getSubmissionStatus(submission)
	if status.status == 'OPEN':

		if externalBucket[0] is not None:
			print 'external bucket %s' % externalBucket[0]
### Change this to give submission object to qsub after get the ability to getSubmissions, downloadFile = False
			cmd = ' '.join(['qsub -N count -o', os.path.join(logsDir, '$JOB_NAME.$JOB_ID'), '-j y -S /usr/bin/python -V', evalCode, '--input', submission.entityId, '--output', outputProjectID[0], '--bucket', externalBucket[0]])
		else:
			cmd = ' '.join(['qsub -N count -o', os.path.join(logsDir, '$JOB_NAME.$JOB_ID'), '-j y -S /usr/bin/python -V', evalCode, '--input', submission.entityId, '--output', outputProjectID])

		print '%s' % cmd
#			subprocess.call(cmd, shell = True)

#		status.status = 'SCORED' # Scored is functioning as "pending" for now.
#		status = syn.store(status)

#if __name__ == "__main__":
	#put test code here?				