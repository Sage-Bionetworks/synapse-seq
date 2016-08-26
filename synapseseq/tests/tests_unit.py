from unittest import TestCase

import synapseseq
import synapseclient
import os

syn = synapseclient.login()


class synio(TestCase):
    def test_command_file(self):
        c = synapseseq.synapse_io.makeCommandFile('funny_prefix',os.getcwd(),syn)
        self.assertTrue(isinstance(c, file))
