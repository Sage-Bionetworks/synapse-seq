"""
Functions for calling quantification methods for RNAseq data.

"""

import os
import subprocess


def runSailfish(prefix,wd,R1,R2,commandsFile=cfHandle):
	'''Runs sailfish on input FASTQ.'''

	if config['workflow']['paired'] is True:	
		cmd = ' '.join(['sailfish quant -i', config['sailfish']['reference'] '-1', R1, '-2', R2,  '-o', os.path.join(wd,prefix), config['sailfish']['options']])
	else:
		cmd = ' '.join(['sailfish quant -i', config['sailfish']['reference'] '-1', R1,  '-o', os.path.join(wd,prefix), config['sailfish']['options']])
	
	expectedResult = os.path.join(wd,prefix,'_'.join([prefix,'quant.sf']))	
	if not os.path.exists(expectedResult):
		print >> commandsFile, '%s' % cmd
		if args.dryrun is False
			subprocess.call(cmd, shell = True)
			os.rename(os.path.join(wd,prefix,'quant.sf'), expectedResult)
	return(expectedResult)
	
