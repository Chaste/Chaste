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

import cppwg.templates.custom


class PopulationWriterCustomTemplate(cppwg.templates.custom.Custom):

    # The writer types this generator instantiates via AddCellWriter<...> etc.
    # These names appear only in the generated code below, never in the parsed
    # C++ signatures, so cppwg cannot auto-include their headers - the generator
    # supplies them itself via get_source_includes().
    cell_writers = [
        "CellAgesWriter",
        "CellAncestorWriter",
        "CellAppliedForceWriter",
        "CellCycleModelProteinConcentrationsWriter",
        "CellDataItemWriter",
        "CellDeltaNotchWriter",
        "CellIdWriter",
        "CellLabelWriter",
        "CellLocationIndexWriter",
        "CellMutationStatesWriter",
        "CellProliferativePhasesWriter",
        "CellProliferativeTypesWriter",
        "CellRadiusWriter",
        "CellRosetteRankWriter",
        "CellVolumesWriter",
        "LegacyCellProliferativeTypesWriter",
    ]

    cell_population_count_writers = [
        "CellMutationStatesCountWriter",
        "CellProliferativeTypesCountWriter",
        "CellProliferativePhasesCountWriter",
    ]

    cell_population_event_writers = [
        "CellDivisionLocationsWriter",
        "CellRemovalLocationsWriter",
    ]

    population_writers = [
        "BoundaryNodeWriter",
        "CellPopulationAdjacencyMatrixWriter",
        "CellPopulationAreaWriter",
        "CellPopulationElementWriter",
        "HeterotypicBoundaryLengthWriter",
        "NodeLocationWriter",
        "NodeVelocityWriter",
        "RadialCellDataDistributionWriter",
        "VertexIntersectionSwapLocationsWriter",
        "VertexT1SwapLocationsWriter",
        "VertexT2SwapLocationsWriter",
        "VertexT3SwapLocationsWriter",
        "VoronoiDataWriter",
    ]

    def __init__(self):
        pass

    def all_writers(self):
        """Return every writer type this generator references."""
        return (
            self.cell_writers
            + self.cell_population_count_writers
            + self.cell_population_event_writers
            + self.population_writers
        )

    def get_source_includes(self, *args, **kwargs):
        """
        Return the writer headers the generated AddCellWriter<...> code needs.

        These are the headers previously listed by hand under the
        population_includes source_includes anchor; keeping them here means the
        writer list and its includes live in one place.
        """
        return [f"{writer}.hpp" for writer in self.all_writers()]

    def get_class_cpp_def_code(self, class_name):
        """
        Custom wrapper code for adding writers to cell populations.
        """

        code = ""

        # Add cell writers
        cell_writer_template = """\
        .def("AddCellWriter{writer}", &{class_name}::AddCellWriter<{writer}>)
"""

        for writer in self.cell_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += cell_writer_template.format(**replacements)

        # Add cell population count writers
        cell_population_count_writer_template = """\
        .def("AddCellPopulationCountWriter{writer}", &{class_name}::AddCellPopulationCountWriter<{writer}>)
"""
        for writer in self.cell_population_count_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += cell_population_count_writer_template.format(**replacements)

        # Add cell population event writers
        cell_population_event_writer_template = """\
        .def("AddCellPopulationEventWriter{writer}", &{class_name}::AddCellPopulationEventWriter<{writer}>)
"""
        for writer in self.cell_population_event_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += cell_population_event_writer_template.format(**replacements)

        # Add population writers
        population_writer_template = """\
        .def("AddPopulationWriter{writer}", &{class_name}::AddPopulationWriter<{writer}>)
"""
        for writer in self.population_writers:
            replacements = {"class_name": class_name, "writer": writer}
            code += population_writer_template.format(**replacements)

        return code
