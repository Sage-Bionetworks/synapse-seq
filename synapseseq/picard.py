"""
Functions for calling Picard.

"""

import os
import subprocess


def samToFastq(localInputFilePath,wd,commandsFile=cfHandle):

	prefix = os.path.basename(localInputFilePath).rstrip('.bam')

	if config['workflow']['paired'] is False:
		R1file = os.path.join(wd, prefix + '_SE.fastq')
		cmd = ' '.join(['java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SortSam.jar'), 'INPUT=', localInputFilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'TMP_DIR=', os.path.join(wd, 'tmp'), 'VALIDATION_STRINGENCY=SILENT'])

		print >> commandsFile, '%s' % cmd
		if args.dryrun is False
			subprocess.call(cmd, shell = True)
		return(R1file)
		
	else:
		R1file = os.path.join(wd, prefix + '_R1.fastq')
		R2file = os.path.join(wd, prefix + '_R2.fastq')
		cmd = ' '.join(['java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SortSam.jar'), 'INPUT=', localInputFilePath, 'OUTPUT=/dev/stdout SORT_ORDER=queryname QUIET=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 TMP_DIR=', os.path.join(wd, 'tmp'), '| java -Xmx2G -jar', os.path.join(config['system']['picard'], 'SamToFastq.jar'), 'INPUT=/dev/stdin FASTQ=', R1file, 'SECOND_END_FASTQ=', R2file, 'TMP_DIR=', os.path.join(wd, 'tmp'), 'VALIDATION_STRINGENCY=SILENT'])

		print >> commandsFile, '%s' % cmd
		if args.dryrun is False
			subprocess.call(cmd, shell = True)
		return((R1file,R2file))	
