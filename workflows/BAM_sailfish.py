#! /usr/bin/env python
# 29 Jul 2016
# Kristen K Dang for Sage Bionetworks


import os
import argparse
import yaml
import synapseclient
import synapseseq
from synapseclient import File, Table, Activity
from toil.job import Job


#################
# Workflow functions
#################

def processInputs(job, args.idFile):

	job.fileStore.logToMaster("Opening list of inputs %s" % args.idFile)
	with open(args.idFile, 'r') as input:
		for line in input:
			job.addChildJobFn(makeFastQ,id=line.strip())
	
	

def makeFastQ(job, id, memory="12G",cores=1,disk="30G"):
	'''Gets input file from Synapse to compute instance.'''

	if args.submission is True:
		submission = synseq.getSubmission(id,args.wd,syn)
	else:
		submission = syn.get(id, downloadLocation=args.wd)
	prefix = submission.name.rstrip('.bam')
	job.fileStore.logToMaster("Downloaded BAM file id: %s, name:" % (submission.id, submission.name))

	cfHandle = synseq.makeCommandFile(prefix,args.wd,syn)
	
	job.fileStore.logToMaster("Running samtoFastQ: %s" % submission.id)
	(R1file, R2file) = samToFastq(os.path.join(args.wd, submission.name),args.wd,commandsFile=cfHandle)
	
	return job.addFollowOnJobFn(quantReads, R1file=R1file, R2file=R2file, submission=submission, cfHandle=cfHandle)


	
	
def quantReads(job, R1file,R2file,prefix,submission,cfHandle, memory="3G",cores=config['workflow']['threads'],disk="500M"):

	job.fileStore.logToMaster("Running sailfish: %s" % submission.id)
	sailResultsFile = runSailfish(prefix=prefix,wd=args.wd,R1=R1file,R2=R2file,commandsFile=cfHandle)
#	return(sailResultsFile)

	## Load output file to Synapse 
	job.fileStore.logToMaster("Loading %s to Synapse" % sailResultsFile)
	sailResultsEntity = File(path=sailResultsFile, name=os.path.basename(sailResultsFile), description='Gene model quant using sailfish.', parentId=config['sailfish']['output'], synapseStore=True)
	sailResultsEntity = syn.store(sailResultsEntity, forceVersion=False)
	syn.setAnnotations(sailResultsEntity, annotations=config['sailfish']['annotations'])
	job.fileStore.logToMaster("New entity id %s" % sailResultsEntity.id)

	## Set up provenance
	job.fileStore.logToMaster("Setting up provenance: %s" % submission.id)
	cfEntity = synseq.storeCommandFile(prefix=prefix,wd=args.wd,syn)

	sailfishAct = Activity(name='Sailfish', description='Quantitate sequencing reads to gene models.', executed=[cfEntity.id])
	sailfishAct.used({'reference':{'targetId':submission.entityId, 'targetVersionNumber':submission.versionNumber}})
	sailfishAct.used(config['sailfish']['reference'])

	syn.setProvenance(sailResultsEntity, sailfishAct) 
	



def changeStatus(submission):

	status = syn.getSubmissionStatus(submission)
	status.status = 'SCORED' # Scored is functioning as "finished" for now.
	status = syn.store(status)



##################
# Initialize
##################


if __name__=="__main__":

	parser = argparse.ArgumentParser(description='Manages sailfish quant pipeline starting with BAM input.')
	parser.add_argument('idFile', help='File containing newline-delimited Synapse id of input files or submissions.')
	parser.add_argument('--conf', required=True, help='Path to local config file.')
	parser.add_argument('--submission', help='Are these submission IDs (not entity)?', action='store_true', default=False)
#	parser.add_argument('--eval', required=False, help='Synapse ID of the eval to run.')
	parser.add_argument('--dryrun', required=False, help='Print all commands that would be executed, but do not run.', action='store_true')
	args = parser.parse_args()


	## Load workflow params
	config = yaml.load(file(args.conf))


	## Log in to Synapse
	syn = synapseclient.login()

		
	## Start workflow	
	Job.Runner.startToil(Job.wrapJobFn(processInputs, args.idFile))