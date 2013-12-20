from distutils.core import setup

setup(name='synapseSeq', version='0.1adev', author='Kristen Dang',
      author_email='kristen.dang@sagebase.org',
      packages=['synapseseq',],
      scripts=['scripts/add_eval.py', 'scripts/eval_browserCoverageTrack.py', 'scripts/eval_fusion_qsub.py', 'scripts/eval_listener.py', 'scripts/eval_raw_counts_htseq.py', 'scripts/eval_SNV_varscan.py'],
      url='https://github.com/Sage-Bionetworks/synapse-seq',
      license='LICENSE.txt',
      description='For loading/running seq data in Synapse.',
      long_description=open('README.md').read(),
      install_requires=['synapseclient>0.5.0'],)
      
