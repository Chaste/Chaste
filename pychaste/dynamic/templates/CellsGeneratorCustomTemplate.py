"""Copyright (c) 2005-2024, University of Oxford.
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


class CellsGeneratorCustomTemplate(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_pre_code(self, class_name):
        """
        Creates custom pybind11 trampoline wrapper code for CellsGenerator.
        """

        code = f"""
class {class_name}_Overrides : {class_name}{{
    public:
    using {class_name}::{class_name};
    std::vector<CellPtr> GenerateBasic(unsigned numCells, const std::vector<unsigned> locationIndices=std::vector<unsigned>(), boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateBasic(cells, numCells, locationIndices, pCellProliferativeType);
        return cells;
    }};
    std::vector<CellPtr> GenerateBasicRandom(unsigned numCells, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateBasicRandom(cells, numCells, pCellProliferativeType);
        return cells;
    }};
    std::vector<CellPtr> GenerateGivenLocationIndices(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateGivenLocationIndices(cells, locationIndices, pCellProliferativeType);
        return cells;
    }};
}};

"""
        return code

    def get_class_cpp_def_code(self, class_name):
        """
        Creates custom pybind11 class wrapper code for CellsGenerator.
        """

        code = f"""
        .def("GenerateBasic",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, const std::vector<unsigned>, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateBasic, " " , py::arg("numCells"),  py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateBasicRandom",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateBasicRandom, " " , py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateGivenLocationIndices",
            (std::vector<CellPtr>({class_name}::*)(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateGivenLocationIndices, " " , py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
"""
        return code
