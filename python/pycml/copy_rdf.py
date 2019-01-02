#!/usr/bin/env python
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

"""
This script copies RDF annotations from a source CellML model to one which has
had all RDF and cmeta:id attributes stripped (e.g. by COR).  It works by
reading the RDF from the original file, and for any variable in the original
with a cmeta:id, adding that id to a variable with the same name and component
in the new model (if any) and adding linked RDF.
"""

import os
import sys

# Make sure PyCml is on sys.path
pycml_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [pycml_path]

import cellml_metadata
import translators

def LoadModel(model_filename, options=[]):
    args = ['-C', '-A', '--assume-valid', model_filename] + options
    options, model_file = translators.get_options(args)
    doc = translators.load_model(model_file, options)
    return doc

def Run(originalModelPath, newModelPath, outputModelPath):
    original_model = LoadModel(originalModelPath)
    new_model = LoadModel(newModelPath)
    # Copy all cmeta:ids to the new model
    vars = original_model.model.xml_xpath(u'cml:component/cml:variable[@cmeta:id]')
    for var in vars:
        try:
            new_var = new_model.model.get_variable_by_name(var.component.name, var.name)
        except:
            continue # No matching variable in the new model
        meta_id = var.cmeta_id
        new_var.xml_set_attribute((u'cmeta:id', translators.NSS['cmeta']), meta_id)
    # Also copy the model's cmeta:id, if any
    model_id = original_model.model.getAttributeNS(translators.NSS['cmeta'], u'id')
    if model_id:
        new_model.model.xml_set_attribute((u'cmeta:id', translators.NSS['cmeta']), model_id)
    # Get the RDF from the original model, and copy to the new one
    for triple in cellml_metadata.get_all_rdf(original_model.model):
        cellml_metadata.add_statement(new_model.model, *triple)
    cellml_metadata.update_serialized_rdf(new_model.model)
    # Write out to a new file
    stream = translators.open_output_stream(outputModelPath)
    new_model.xml(indent=u'yes', stream=stream)
    translators.close_output_stream(stream)
    # Tidy up nicely
    cellml_metadata.remove_model(original_model.model)
    cellml_metadata.remove_model(new_model.model)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print >> sys.stderr, "Usage:", sys.argv[0], "<original_model> <new_model> <output_path>"
        sys.exit(1)
    Run(*sys.argv[1:])
