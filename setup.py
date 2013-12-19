from distutils.core import setup

setup(name='synapse-seq', version='0.1a', author='Kristen Dang',
      author_email='kristen.dang@sagebase.org',
      packages='synapse-seq',
      scripts=['scripts/add_eval.py, scripts/eval_browserCoverageTrack.py, scripts/eval_fusion_qsub.py, scipts/eval_listener.py, scripts/eval_raw_counts.htseq.py, scripts/eval_SNV_varscan.py'],
      url='https://github.com/Sage-Bionetworks/synapse-seq',
      license='LICENSE.txt',
      description='For loading/running seq data in Synapse.',
      long_description=open('README.md')read(),
      install_requires=['synapseclient>0.5.0'],)
      
