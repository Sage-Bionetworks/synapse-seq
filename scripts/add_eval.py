#! /usr/bin/env python
# KKD for Sage Bionetworks
# Create synapse evaluation entity

import synapseclient, argparse
from synapseclient import Evaluation

parser = argparse.ArgumentParser(description='Creates XYZ evaluation via synapse.')
parser.add_argument('--name', dest='name', required=True, help='Name of workflow')
parser.add_argument('--description', dest='desc', required=True, help='Description of workflow')
parser.add_argument('--parent', dest='pid', required=False, help='Project ID with which to associate this evaluation.', default='syn1972151')

args = parser.parse_args()


syn = synapseclient.Synapse()
syn.login()

newEval = Evaluation(name = args.name, description = args.desc, contentSource = args.pid)
newEval = syn.store(newEval)
syn.joinEvaluation(newEval)

print 'Synapse ID for evaluation %s is %s' % (args.name, newEval.id)
