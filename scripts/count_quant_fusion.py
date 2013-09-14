#! /usr/bin/env python
# KKD for Sage Bionetworks
# Aug. 16, 2013


import argparse, os, subprocess, sys

parser = argparse.ArgumentParser(description='Count-quant-fusion pipeline for seq workflows.')
parser.add_argument('--bamdir', dest='bam', required=True, help='Location of BAM files to process.')

args = parser.parse_args()



# Write pipeline steps to file and qsub
file_list = os.listdir(args.bam) 
for file in file_list:
	
	if file.endswith('bam'):
		stepIDs = dict()

		jobFileName = ''.join(['tempfile_', file.rstrip('.bam'), '.sh'])
		jobfile = open(jobFileName, 'w')

		# write header		
		jobfile.write('#! /bin/bash\n\n')
		jobfile.write('#$ -N cqf\n#$ -V\n')
		jobfile.write('#$ -o /home/kkdang/logs/$JOB_NAME.$JOB_ID -j y\n\n')

		# load programs
		jobfile.write('module load samtools picard-tools bowtie cufflinks\n\n')
			
		# write samsort step
		filePrefix = file.rstrip('.bam') + '_namesort'				
		jobfile.write(' '.join(['samtools sort -n', file, filePrefix, '\n']))
			
		# write htseq step
		sortedFile = filePrefix + '.bam'				
		outputFile = sortedFile.rstrip('.bam') + '.htseq'				
		jobfile.write(' '.join(['samtools view ', sortedFile, '| python -m HTSeq.scripts.count -m intersection-nonempty -s no - /home/kkdang/data/ccle/Hsapiens_Ensembl_v70_noPNH.gtf >', outputFile, '\n']))

		# write dexseq step		
		outputFile = sortedFile.rstrip('.bam') + '.dexseq'				
		jobfile.write(' '.join(['samtools view', sortedFile, '| python ~/bin/dexseq_count.py -p yes -s no /external-data/Genome/gene_models/Hsapiens_Ensembl_v70_dex.gtf -', outputFile, '\n']))
	
		# write cufflinks no-assembly step
		filePrefix = file.rstrip('.bam') + '_namesort'	
		outdir = os.path.join(args.bam, file.rstrip('.bam'))
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		jobfile.write(' '.join(['cufflinks --output-dir', outdir, '--GTF /home/kkdang/data/ccle/Hsapiens_Ensembl_v70_noPNH.gtf --frag-bias-correct /external-data/Genome/genomes/Hsapiens_Ensembl_GRCh37/WholeGenomeFasta/genome.fa --multi-read-correct --compatible-hits-norm', file, '\n\n']))

		# write samtofastq step
		## *** need to finish this command
		R1file = file.rstrip('.bam') + '_R1.fastq'
		R2file = file.rstrip('.bam') + '_R2.fastq'
		outputFile = file.rstrip('.bam') + '_trans.bam'
		jobfile.write(' '.join(['/gluster/toolbox/picard-tools/1.94/SamToFastq.jar\n']))

		# write bowtie2-eXpress step
		jobfile.write(' '.join(['bowtie2 --time -p 4 -k 1000 --very-sensitive -x Hsapiens_Ensembl_v70 -1', R1file, '-2', R2file, '| /home/kkdang/bin/expr4.0/express /external-data/Genome/indicies/Hsapiens_Ensembl_v70GTF_Bowtie2/Hsapiens_Ensembl_v70.fa\n']))

		# write TIGAR step
		## *** need to finish this command
		jobfile.write(' '.join(['bowtie2 --time -p 4 -k 1000 --very-sensitive -x Hsapiens_Ensembl_v70 -1', R1file, '-2', R2file, '| samtools view -Sb - >', outputFile, '\n']))
		filePrefix = file.rstrip('.bam') + '_namesort'				
		jobfile.write(' '.join(['TIGAR', file, filePrefix]))



		# submit job
		jobfile.close()
		cmd = " ".join(['qsub', jobFileName])
		print '%s' % cmd
	#	submit = call("qsub submit_dexseq.sh jobFileName", shell=True)

