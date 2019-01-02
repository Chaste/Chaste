
"""Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
import unittest

# Get PyCml modules
import cellml_metadata
import translators

def LoadModel(model_filename, options=[]):
    args = ['-C', '-A', '--assume-valid', model_filename] + options
    options, model_file = translators.get_options(args)
    doc = translators.load_model(model_file, options)
    return doc

class TestCopyingRdf(unittest.TestCase):
    def _CompareRdf(self, referenceFilePath, outputFilePath):
        """Helper method to compare the RDF content of two CellML files."""
        ref_model = LoadModel(referenceFilePath)
        output_model = LoadModel(outputFilePath)
        ref_rdf = list(cellml_metadata.get_all_rdf(ref_model.model))
        output_rdf = list(cellml_metadata.get_all_rdf(output_model.model))
        self.assertEqual(len(ref_rdf), len(output_rdf))
        # Tidy up for Redland
        cellml_metadata.remove_model(ref_model.model)
        cellml_metadata.remove_model(output_model.model)
    
    def TestCopyingRdfOnLr91(self):
        original_model = 'heart/src/odes/cellml/LuoRudy1991.cellml'
        no_rdf_model = 'python/test/data/LuoRudy1991WithNoRdf.cellml'
        output_model = os.path.join(CHASTE_TEST_OUTPUT, 'LuoRudy1991Output.cellml')
        rc = os.system('python/pycml/copy_rdf.py ' + original_model + ' ' + no_rdf_model + ' ' + output_model)
        self.assertEqual(rc, 0)
        # Compare output to reference
        self._CompareRdf('python/test/data/LuoRudy1991Output.cellml', output_model)
