/*

Copyright (c) 2005-2024, University of Oxford.
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

*/

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"

#include "CellCycleModelProteinConcentrationsWriter_3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellCycleModelProteinConcentrationsWriter<3, 3> CellCycleModelProteinConcentrationsWriter_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellCycleModelProteinConcentrationsWriter_3_3_Overrides : public CellCycleModelProteinConcentrationsWriter_3_3
{
public:
    using CellCycleModelProteinConcentrationsWriter_3_3::CellCycleModelProteinConcentrationsWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            double,
            CellCycleModelProteinConcentrationsWriter_3_3,
            GetCellDataForVtkOutput,
            pCell,
            pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            CellCycleModelProteinConcentrationsWriter_3_3,
            VisitCell,
            pCell,
            pCellPopulation);
    }
};

void register_CellCycleModelProteinConcentrationsWriter_3_3_class(py::module &m)
{
    py::class_<CellCycleModelProteinConcentrationsWriter_3_3, CellCycleModelProteinConcentrationsWriter_3_3_Overrides, boost::shared_ptr<CellCycleModelProteinConcentrationsWriter_3_3>, AbstractCellWriter<3, 3>>(m, "CellCycleModelProteinConcentrationsWriter_3_3")
        .def(py::init<>())
        .def("GetCellDataForVtkOutput",
            (double(CellCycleModelProteinConcentrationsWriter_3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellCycleModelProteinConcentrationsWriter_3_3::GetCellDataForVtkOutput,
            " ", py::arg("pCell"), py::arg("pCellPopulation"))
        .def("VisitCell",
            (void(CellCycleModelProteinConcentrationsWriter_3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellCycleModelProteinConcentrationsWriter_3_3::VisitCell,
            " ", py::arg("pCell"), py::arg("pCellPopulation"))
    ;
}
