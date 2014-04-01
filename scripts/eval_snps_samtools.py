#! /usr/bin/env python
# Kristen K. Dang for Sage Bionetworks
# Mar. 21, 2014


import synapseclient, os, argparse, subprocess, math, shutil, sys, hashlib
from synapseclient import File

parser = argparse.ArgumentParser(description='Runs samtools variant calling on RNAseq data.')
parser.add_argument('--submission', required=True, dest='bam', help='Submission id of BAM file to process.')
parser.add_argument('--params', required=True, help='Synapse ID of the parameter file for this job.')
parser.add_argument('--ref', required=True, help='Synapse ID or public URL of reference genome used to do the alignment.', default=None)
parser.add_argument('--output', dest='out', required=False, help='Synapse ID of project/folder to which results are uploaded.', default='syn2384001')
parser.add_argument('--limitPositions', dest='limit', required=False, help='Path to file of positions at which to call variants, if not doing genome-wide.', default=None)
parser.add_argument('--depth', required=False, help='Constant by which average read depth is multiplied to get the maximum depth.', default=50)
parser.add_argument('--bucket', dest='bucket', required=False, help='Bucket name if not stored in Synapse.', default=None)
parser.add_argument('--keyname', dest='key', required=False, help='Key name if not stored in Synapse.', default=None)
args = parser.parse_args()


# Hardcoded to cloudbiolinux - need to generalize eventually
bcftools = '/usr/bin/bcftools' # assume location of cloudbiolinux
varFilter = '/usr/share/samtools/vcfutils.pl varFilter' # assume location of cloudbiolinux
reference = '/mnt/transient_nfs/genome.fa' # assume location of cloudbiolinux master node
if not os.path.exists('home/ubuntu/.synapseConfig'):
	shutil.copy('/mnt/transient_nfs/.synapseConfig', '/home/ubuntu/')
wd = '/mnt/galaxyData/samtools_results'
if not os.path.exists(wd):
	subprocess.call('sudo chown ubuntu:ubuntu /mnt/galaxyData', shell = True)
	os.mkdir(wd)

def calc_md5(inFile): # This function code from Chris Bare
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
	
	
syn = synapseclient.Synapse()
syn.login()
submission = syn.getSubmission(args.bam, downloadFile = False)

if args.bucket is not None:
	import boto
	s3 = boto.connect_s3()
	external_bucket = s3.get_bucket(args.bucket) 
	bucketItem = external_bucket.get_key(args.key)
	localBAMfilePath = os.path.join(wd, os.path.basename(args.key))
	if not os.path.exists(localBAMfilePath):
		print 'Getting data in bucket %s' % args.bucket
		bucketItem.get_contents_to_filename(localBAMfilePath)	
else:
	localBAMfilePath = os.path.join(wd, submission.name)
	if not os.path.exists(localBAMfilePath):
		print 'Will download %s from synapse.' % submission.name
		BAMentity = syn.get(entity=submission.entityId, version=submission.versionNumber, downloadFile = True, downloadLocation = wd)
		if calc_md5(localBAMfilePath).upper() != BAMentity.md5.upper():
			print 'Download error for file %s' % submission.name
			os.remove(localBAMfilePath)
			sys.exit()


print '%s' % localBAMfilePath	
prefix = os.path.basename(localBAMfilePath).rstrip('.bam')

## Run variant calling step
outRawBCF = os.path.join(wd, prefix+'.var.raw.bcf')
if args.limit is not None:
	cmd = ' '.join(['samtools mpileup -uf', reference, '-l', args.limit, localBAMfilePath, '|  bcftools view -bvcg - >', outRawBCF])
else:
	cmd = ' '.join(['samtools mpileup -uf', reference, localBAMfilePath, '|  bcftools view -bvcg - >', outRawBCF])

print '%s' % cmd
subprocess.call(cmd, shell = True)


## Run variant filtering step
# calculate the average read depth for this sample
cmd  = ' '.join([bcftools, 'view', outRawBCF, '| awk \'$1 !~ /^#/ {print $8}\' | awk --field-separator=\';\' \'$1 ~ /DP=/ {print $1}\' | cut -f2 -d= | awk \'{sum+=$1}END{print sum/NR}\''])
depth = subprocess.check_output(cmd, shell = True).split()[0]

outFltVCF = os.path.join(wd, prefix+'.var.flt.vcf')
cmd = ' '.join(['bcftools view', outRawBCF,  '|', varFilter, '-D'+str(math.ceil(args.depth*float(depth))), '>', outFltVCF])
print '%s' % cmd
subprocess.call(cmd, shell = True)


## Load filtered results to synapse
print 'Loading %s to Synapse.' % outFltVCF

vcf = File(path=outFltVCF, name=os.path.basename(outFltVCF), description='Filtered variant calls.', parentId=args.out, synapseStore=True)	

vcf = syn.store(vcf, activityName='variant calling', activityDescription='Default variant calling on RNAseq data.', forceVersion=False, used=[BAMentity.id, args.ref, args.params], executed = ['syn2243148', 'https://github.com/Sage-Bionetworks/synapse-seq/blob/master/scripts/eval_snps_samtools.py'])

syn.setAnnotations(vcf, annotations=dict(fileType='VCF',VariantOnly='True',LimitedBy=args.limit))
print 'new entity id %s' % vcf.id


# clean up BAM files
os.remove(localBAMfilePath)

# change status of BAM 
status = syn.getSubmissionStatus(submission)
status.status = 'SCORED' # Scored is functioning as "pending" for now.
status = syn.store(status)
