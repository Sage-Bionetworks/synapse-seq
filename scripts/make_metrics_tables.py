#! /usr/bin/env python
# June 30, 2015
# Kristen K Dang for Sage Bionetworks


import os, argparse
import synapseclient
from synapseclient import Schema, Column

parser = argparse.ArgumentParser(description='Creates empty metrics tables for a seq workflow.')
parser.add_argument('synProject', help='Synapse id of project that will contain the table(s).')
parser.add_argument('--star', required=False, dest='star', help='Create table for metrics output from star.', action='store_true')
parser.add_argument('--featurecounts', required=False, dest='fc', help='Create table for metrics output from featurecounts.', action='store_true')
parser.add_argument('--sailfish', required=False, dest='sailfish', help='Create table for metrics output from sailfish.', action='store_true')
parser.add_argument('--picardRNA', required=False, dest='picard', help='Create table for metrics output from picardRNA.', action='store_true')
parser.add_argument('--trimming', required=False, dest='trim', help='Create table for metrics output from fastq trimming.', action='store_true')
args = parser.parse_args()

syn = synapseclient.login()

if args.star is True:

	cols = [
		Column(name='Sample', columnType='STRING', maximumSize=40),
		Column(name='SubmissionID', columnType='STRING', maximumSize=12),
		Column(name='Started-job', columnType='DATE'),
		Column(name='Started-map', columnType='DATE'),
		Column(name='Finished', columnType='DATE'),
		Column(name='Mapping-speed', columnType='DOUBLE'),
		Column(name='Number-reads', columnType='INTEGER'),
		Column(name='Ave-read-length', columnType='DOUBLE'),
		Column(name='Uniq-reads', columnType='INTEGER'),
		Column(name='Uniq-reads-PCT', columnType='DOUBLE'),
		Column(name='Ave-mapped-len', columnType='DOUBLE'),
		Column(name='Total-splices', columnType='INTEGER'),
		Column(name='Annotated-splices', columnType='INTEGER'),
		Column(name='GT/AG-splices', columnType='INTEGER'),
		Column(name='GC/AG-splices', columnType='INTEGER'),
		Column(name='AT/AC-splices', columnType='INTEGER'),
		Column(name='Non-canonical-splices', columnType='INTEGER'),
		Column(name='Mismatch-rate-per-base-PCT', columnType='DOUBLE'),
		Column(name='Deletion-rate-per-base', columnType='DOUBLE'),
		Column(name='Deletion-ave-length', columnType='DOUBLE'),
		Column(name='Insertion-rate-per-base', columnType='DOUBLE'),
		Column(name='Insertion-ave-length', columnType='DOUBLE'),
		Column(name='Reads-mapped-to-multiple-loci', columnType='INTEGER'),
		Column(name='Reads-mapped-to-multiple-loci-PCT', columnType='DOUBLE'),
		Column(name='Reads-mapped-to-too-many-loci', columnType='INTEGER'),
		Column(name='Reads-mapped-to-too-many-loci-PCT', columnType='DOUBLE'),
		Column(name='Reads-unmapped-too-many-mismatches-PCT', columnType='DOUBLE'),
		Column(name='Reads-unmapped-too-short-PCT', columnType='DOUBLE'),
		Column(name='Reads-unmapped-other-PCT', columnType='DOUBLE')
	]

	schema = syn.store(Schema(name='STAR_alignment_metrics', columns=cols, parent=args.synProject))
	

if args.picard is True:

	cols = [
		Column(name='Sample', columnType='STRING', maximumSize=40),
		Column(name='SubmissionID', columnType='STRING', maximumSize=12),
		Column(name='PF_BASES', columnType='INTEGER'),
		Column(name='PF_ALIGNED_BASES', columnType='INTEGER'),
		Column(name='RIBOSOMAL_BASES', columnType='INTEGER'),
		Column(name='CODING_BASES', columnType='INTEGER'),
		Column(name='UTR_BASES', columnType='INTEGER'),
		Column(name='INTRONIC_BASES', columnType='INTEGER'),
		Column(name='INTERGENIC_BASES', columnType='INTEGER'),
		Column(name='IGNORED_READS', columnType='INTEGER'),
		Column(name='CORRECT_STRAND_READS', columnType='INTEGER'),
		Column(name='INCORRECT_STRAND_READS', columnType='INTEGER'),
		Column(name='PCT_RIBOSOMAL_BASES', columnType='DOUBLE'),
		Column(name='PCT_CODING_BASES', columnType='DOUBLE'),
		Column(name='PCT_UTR_BASES', columnType='DOUBLE'),
		Column(name='PCT_INTRONIC_BASES', columnType='DOUBLE'),
		Column(name='PCT_INTERGENIC_BASES', columnType='DOUBLE'),
		Column(name='PCT_MRNA_BASES', columnType='DOUBLE'),
		Column(name='PCT_USABLE_BASES', columnType='DOUBLE'),
		Column(name='PCT_CORRECT_STRAND_READS', columnType='DOUBLE'),
		Column(name='MEDIAN_CV_COVERAGE', columnType='DOUBLE'),
		Column(name='MEDIAN_5PRIME_BIAS', columnType='DOUBLE'),
		Column(name='MEDIAN_3PRIME_BIAS', columnType='DOUBLE'),
		Column(name='MEDIAN_5PRIME_TO_3PRIME_BIAS', columnType='DOUBLE'),
		Column(name='SAMPLE', columnType='STRING', maximumSize=40),
		Column(name='LIBRARY', columnType='STRING', maximumSize=40),
		Column(name='READ_GROUP', columnType='STRING', maximumSize=40)
	]

	schema = syn.store(Schema(name='Picard_RNA_alignment_metrics', columns=cols, parent=args.synProject))



if args.fc is True:

	cols = [
		Column(name='Sample', columnType='STRING', maximumSize=40),
		Column(name='SubmissionID', columnType='STRING', maximumSize=12),
		Column(name='Assigned', columnType='INTEGER'),
		Column(name='Unassigned_Ambiguity', columnType='INTEGER'),
		Column(name='Unassigned_MultiMapping', columnType='INTEGER'),
		Column(name='Unassigned_NoFeatures', columnType='INTEGER'),
		Column(name='Unassigned_Unmapped', columnType='INTEGER'),
		Column(name='Unassigned_MappingQuality', columnType='INTEGER'),
		Column(name='Unassigned_FragmentLength', columnType='INTEGER'),
		Column(name='Unassigned_Chimera', columnType='INTEGER'),
		Column(name='Unassigned_Secondary', columnType='INTEGER')
	]

	schema = syn.store(Schema(name='featurecounts_metrics', columns=cols, parent=args.synProject))


if args.trim is True:

	cols = [
		Column(name='Sample', columnType='STRING', maximumSize=40),
		Column(name='SubmissionID', columnType='STRING', maximumSize=12),
		Column(name='R1_file', columnType='STRING', maximumSize=70),
		Column(name='R1_num_reads', columnType='INTEGER'),
		Column(name='R1_num_empty_reads', columnType='INTEGER'),
		Column(name='R1_num_trimmed_reads', columnType='INTEGER'),
		Column(name='R1_num_trimmed_bases', columnType='DOUBLE'),
		Column(name='R2_file', columnType='STRING', maximumSize=70),
		Column(name='R2_num_reads', columnType='INTEGER'),
		Column(name='R2_num_empty_reads', columnType='INTEGER'),
		Column(name='R2_num_trimmed_reads', columnType='INTEGER'),
		Column(name='R2_num_trimmed_bases', columnType='DOUBLE')
	]

	schema = syn.store(Schema(name='FASTQ_trim_metrics', columns=cols, parent=args.synProject))
	
