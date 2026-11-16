"""Copyright (c) 2005-2026, University of Oxford.
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

import importlib.util
import os
import sys


def _load_population_writer_template():
    """
    Load PopulationWriterCustomTemplate from this directory.

    cppwg loads each custom generator by file path rather than by import, so
    the templates directory is not on sys.path. Load the base generator the
    same way, under the same module name cppwg uses, so both share one module.
    """
    path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "PopulationWriterCustomTemplate.py",
    )
    name = os.path.splitext(path)[0]
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name].PopulationWriterCustomTemplate


class OffLatticePopulationWriterCustomTemplate(_load_population_writer_template()):

    # AddNodePointDataWriter<...> is declared on AbstractOffLatticeCellPopulation,
    # so this generator is only used for off-lattice populations.
    node_point_data_writers = [
        "BoundaryNodePointDataWriter",
        "ElementIdNodePointDataWriter",
        "NodeRegionPointDataWriter",
    ]

    def all_writers(self):
        """Return every writer type this generator references."""
        return super().all_writers() + self.node_point_data_writers

    def get_class_cpp_def_code(self, class_name):
        """
        Custom wrapper code for adding writers to off-lattice cell populations.
        """

        code = super().get_class_cpp_def_code(class_name)

        # Add node point data writers
        node_point_data_writer_template = """\
        .def("AddNodePointDataWriter_{writer}", &{class_name}::AddNodePointDataWriter<{writer}>)
"""
        for writer in self.node_point_data_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += node_point_data_writer_template.format(**replacements)

        return code
