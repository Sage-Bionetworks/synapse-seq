"""
Functions for Synapse interactions of seq workflows.

"""


import synapseclient
import hashlib
import os


def getBAMs(projectOrFolderID, syn):
	'''Returns list of BAM entities in the given project or folder id.'''
	
	# Check if projectOrFolderID is BAM folder
	container = syn.get(projectOrFolderID, downloadFile = False)
	if container.name == 'BAM':
		BAMfolderID = container.id	
	# If not, get ID of BAM folder for this container
	
	BAMEntityList = list()
	print 'BAM folder ID: %s' % BAMfolderID
	sql = ''.join(['SELECT id, name FROM entity WHERE parentId=="', BAMfolderID, '"'])
	print '%s' % sql
	results = syn.chunkedQuery(sql)
	for result in results:
		if result['entity.name'].endswith('.bam'):
			BAMEntityList.append(syn.get(result['entity.id'], downloadFile = False))
	
	return(BAMEntityList)
	


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



def getSubmittedBAM(evalID,syn): 
	'''Get dictionary of all BAMs submitted to given eval.'''

	submittedEntities = dict()
	for submission in syn.getSubmissions(evalID):	
		BAMentity = syn.get(submission['entityId'], version = submission['versionNumber'], downloadFile = False)	
		submittedEntities[BAMentity.name] = BAMentity.id
#		print '%s' % BAMentity.name
	return(submittedEntities)
	

def getSubmission(submissionID,wd,syn):
	'''Downloads submission entity.'''

	submission = syn.getSubmission(submissionID, downloadFile = False)
	localInputFilePath = os.path.join(wd, submission.name)
	if not os.path.exists(localInputFilePath):
		submission = syn.getSubmission(args.submission, downloadFile = True, downloadLocation = wd)
	return(submission)


def makeCommandFile(prefix,wd,syn):

	cfPath = os.path.join(wd, '_'.join([prefix, 'commands.txt']))
	commandsFile = open(os.path.join(wd, cfPath),'w')
	return(commandsFile)
	
def storeCommandFile(cfHandle,wd,syn):

	cfHandle.close()
	cfEntity = File(path=os.path.join(wd, cfHandle.name), description='Job commands.', parentId=config['workflow']['output'], synapseStore=True)
	cfEntity = syn.store(cfEntity)
	return(cfEntity)

	




#if __name__ == "__main__":
	#put test code here