"""
Functions for loading existing sequencing datasets into Synapse.

"""

# To do items: 
#	change config to yaml
#	add dependencies in provenance commands
#	generalize tophat name/version id for other methods
#	unit tests

from synapseclient import Project, Folder, File, Activity, Evaluation
import subprocess, os


def readConfig(configFile):
	'''Reads external config file and stores entries in dictionary allArgDict organized by workflow step.'''
	
	allArgDict = dict()
	foldersToCreate = list()
	workflow = False
	stepNum = 0
	params = open(configFile, 'r')
	finalStep = 1
	for line in params:
		if line.startswith('# workflow'):
			workflow = True
			if stepNum > 0:
				allArgDict[stepNum] = workflowDict
			stepNum = line.strip().split()[3]
			workflowDict = dict()
			finalStep = finalStep + 1
			continue
		elif line.startswith('#') or line.startswith('\n'): 
			continue
		fields = line.strip().split('\t')
	
		if workflow == True:
			if fields[0] == 'annotations':
				annotDict = dict()
				elements = fields[1].split(';')
				for item in elements:
					key, value = item.split(',')
					annotDict[key] = value
				workflowDict['annotations'] = annotDict
			elif fields[0] == 'folder':
				foldersToCreate.append(fields[1])
				workflowDict[fields[0]] = fields[1]
			else:
				workflowDict[fields[0]] = fields[1]
		else:
			try:
				allArgDict[fields[0]] = fields[1]
			except IndexError: 
				continue
	allArgDict[stepNum] = workflowDict
	allArgDict['folderlist'] = foldersToCreate
	allArgDict['finalStep'] = finalStep
	return(allArgDict)


## test code to make sure dictionary is formed correctly
# for item in allArgDict:
# 	try: 
# 		for k, v in allArgDict[item].iteritems():
# 			print '%s %s' % (k, v)
# 	except AttributeError:
# 		print '%s - not a dictionary' % item 
# 		continue



def parseTophatCall(inBAMPath): 
	'''For BAM files generated with TopHat, parses the cmd statement and tophat version from the sam header.'''
	
	cmd = ' '.join(['samtools view -H', inBAMPath])
	header = subprocess.check_output(cmd, shell=True)
	headerList = header.split('\n')
	for line in headerList:
		if line.startswith('@PG'):
			elements = line.split('\t')
			name = elements[1].split(':')[1]
			version = elements[2].split(':')[1]
#			print '%s %s' % (name, version)
			return(' '.join([name,version]), elements[3])


def str2bool(inString): # this function found on stackoverflow
	return inString.lower() in ("yes", "true", "t", "1")


def add_workflow_step_to_synapse(inFilePath, stepDict, step='1', store=False, software=None, parentid=None, syn=None, stepIDs=None, inFilename=None ):
	'''Uploads files with provenance and annotations to Synapse.'''
	if not inFilename:
		inFilename = os.path.basename(inFilePath.strip())
	if not software:
		software = stepDict['softwareName']
	if 'used' in stepDict:
		usedList = stepDict['used'].strip().split(',')
		if 'depends' in stepDict:
			usedList.append(stepIDs[stepDict['depends']])
	elif 'depends' in stepDict:
		usedList = stepIDs[stepDict['depends']]
	execList = stepDict['executed'].strip().split(';')

	act = Activity(name=stepDict['actName'],  description=stepDict['description'])
	act.used(usedList)
	for item in execList:
		splitItem = item.split(',')
		target = splitItem[0]
		version = 1
		if (len(splitItem) > 1):
			version = splitItem[1]
		act.executed(target=target, targetVersion=version)

	step_file = File(inFilePath, name=inFilename, description=stepDict['fileDescription'], parentId=parentid, synapseStore=str2bool(stepDict['store']))	
	step_file = syn.store(step_file, activity=act)
	syn.setAnnotations(step_file, annotations=stepDict['annotations'])
	return(step_file.id)


def setUpSynapseProject(foldersToCreate,syn,pid=None,pname=None):
	'''Creates Synapse project and necessary folders for the dataset.'''
	
	# Create a set of sub-folders expected for this project
	folderSchemaSet = set(foldersToCreate)
	existingFolders = dict()

	# Get the project if it exists or create
	if pid == None:
		project = Project(pname)
		project = syn.store(project)

	else:
		project = syn.get(pid)
		print '%s' % project.name
		# Check to see whether any of the sub-folders already exist.
		results = syn.chunkedQuery(''.join(['SELECT id, name FROM entity WHERE parentId=="', project.id, '"']))
		# Add the ones that exist
		for entity in results:
			existingFolders[entity['entity.name']] = entity['entity.id']
		# Make a list of ones to create
		if len(existingFolders) > 0:
			foldersToCreate = folderSchemaSet.difference(existingFolders.keys())

	# create the folders that don't exist
	for name in foldersToCreate:
		createFolder = Folder(name, parent=project.id)
		createFolder = syn.store(createFolder)
		existingFolders[name] = createFolder.id
	return(project, existingFolders)


