from distutils.core import setup

setup(name='synapseseq', version='0.2a', author='Kristen Dang',
      author_email='kristen.dang@sagebase.org',
      packages=['synapseseq',],
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['scripts/add_eval.py'],
      url='https://github.com/Sage-Bionetworks/synapse-seq',
      license='LICENSE.txt',
      description='For running seq pipelines on [cloud|local] clusters using Synapse to track data and provenance.',
      long_description=open('README.md').read(),
      install_requires=['synapseclient>1.0','pyyaml'],)
      
