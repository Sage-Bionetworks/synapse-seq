"""
Functions for running sequencing workflows using Synapse.

"""

# TODO unit tests

import synapseclient, os, argparse, sys, subprocess, hashlib, shutil
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



# This function code from Chris Bare
def calc_md5(inFile): 
	'''Calculate md5.'''
	
	md5 = hashlib.md5()
	with open(inFile, 'rb') as bamfile:
		while True:
			data = bamfile.read(2**20) 
			if not data:
				break
			md5.update(data)
	bamfile.closed
	print '%s' % md5.hexdigest().upper()
	return(md5.hexdigest())


def locateRefOnHeadNode(ref,headNFSPath):
	'''Locate reference on head node NFS.'''

	if ref.startswith('syn'):
		refEntity = syn.get(entity=ref, downloadFile = False)
		refName = refEntity.name
	else: # for external links
		refName = os.path.basename(ref)
	reference = os.path.join(headNFSPath, refName)
	if not os.path.exists(reference):
		sys.exit("ERROR: Can't locate reference.")
	return(reference)


def copyRefToWorkerNode(ref,headNFSPath,localPath):
	'''Copy reference from head node NFS to worker.'''
	
	refOnHead = locateRefOnHeadNode(ref=ref,headNFSPath=headNFSPath)
	localRefPath = os.path.join(localPath,os.path.basename(refOnHead))
	if not os.path.exists(localRefPath):
		shutil.copy(refOnHead, localPath)
	return(localRefPath)


def getBAMtoComputeNode(wd,submission=None,bucket=None,extKey=None):
	'''Locate BAM on node, copying from synapse or S3 if necessary.'''

	if bucket is not None:
		import boto
		s3 = boto.connect_s3()
		external_bucket = s3.get_bucket(bucket) 
		bucketItem = external_bucket.get_key(extKey)
		localBAMfilePath = os.path.join(wd, os.path.basename(extKey))
		if not os.path.exists(localBAMfilePath):
			print 'Getting data in bucket %s' % bucket
			bucketItem.get_contents_to_filename(localBAMfilePath)	
	else:
		localBAMfilePath = os.path.join(wd, submission.name)
		if not os.path.exists(localBAMfilePath):
			print 'Will download %s from synapse.' % submission.name
			BAMentity = syn.get(entity=submission.entityId, version=submission.versionNumber, downloadFile = True, downloadLocation = wd)
			if calc_md5(localBAMfilePath).upper() != BAMentity.md5.upper():
				os.remove(localBAMfilePath)
				sys.exit(' '.join(['ERROR: MD5 of %s does not match md5 in synapse.', submission.name]))
		else:
			BAMentity = syn.get(entity=submission.entityId, version=submission.versionNumber, downloadFile = False, downloadLocation = wd)
	return(localBAMfilePath)



def getSynConfigFromHead(localPath,headPath):
	'''Copy Synapse config file from head to worker node, if necessary.'''

	if not os.path.exists(os.path.join(localPath, '.synapseConfig')):
		shutil.copy(os.path.join(headPath, '.synapseConfig'), localPath)


#if __name__ == "__main__":
	#put test code here?				