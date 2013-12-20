synapseseq
===========

This package contains functionality to load existing seq datasets into [Synapse](https://www.synapse.org/) and to run seq workflows using [Synapse.](https://www.synapse.org/)

This software is in early alpha stage -- please contact the author before attempting to use.


Installing
------------------

Clone the git repo, or download the zip, or to install local git repo:
    sudo pip install -e /path/to/local/repo


## Dependencies
synapse python client, qsub


Usage
------------------

### Running seq workflows using synapse

Quick overview:

1. Install synapse python client and synapseseq.
2. Using synapse client, create evaluations for all workflow components you want to use (components are found in scripts/).
3. Edit the file eval-code-assignment.yaml to include the ids for the evals you created in step #2. Upload the file to synapse and edit scripts/eval\_listener to point to this entity.
4. Setup scripts/eval\_listener.py to run on a cron job or other comparable method. 
5. Upload input files to synapse or external AWS S3 bucket and submit them to desired seq workflows (aka evals). Outputs will be loaded to synapse with provenance and annotations.


Compute requirements:

Currently, qsub is used to manage all jobs across the compute resource. Other job management software may be added in the future.


### Loading existing datasets into Synapse.

Will add this text in near future.


License and Copyright
---------------------

&copy; Copyright 2013 Sage Bionetworks

This software is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
