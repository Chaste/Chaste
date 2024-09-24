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
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "DifferentialAdhesionPottsUpdateRule_3.cppwg.hpp"

namespace py = pybind11;
typedef DifferentialAdhesionPottsUpdateRule<3> DifferentialAdhesionPottsUpdateRule_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DifferentialAdhesionPottsUpdateRule_3_Overrides : public DifferentialAdhesionPottsUpdateRule_3
{
public:
    using DifferentialAdhesionPottsUpdateRule_3::DifferentialAdhesionPottsUpdateRule;
    double GetCellCellAdhesionEnergy(::CellPtr pCellA, ::CellPtr pCellB) override
    {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule_3,
            GetCellCellAdhesionEnergy,
            pCellA,
            pCellB);
    }
    double GetCellBoundaryAdhesionEnergy(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule_3,
            GetCellBoundaryAdhesionEnergy,
            pCell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            DifferentialAdhesionPottsUpdateRule_3,
            OutputUpdateRuleParameters,
            rParamsFile);
    }
};

void register_DifferentialAdhesionPottsUpdateRule_3_class(py::module &m)
{
    py::class_<DifferentialAdhesionPottsUpdateRule_3, DifferentialAdhesionPottsUpdateRule_3_Overrides, boost::shared_ptr<DifferentialAdhesionPottsUpdateRule_3>, AdhesionPottsUpdateRule<3>>(m, "DifferentialAdhesionPottsUpdateRule_3")
        .def(py::init<>())
        .def("GetCellCellAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule_3::*)(::CellPtr, ::CellPtr)) &DifferentialAdhesionPottsUpdateRule_3::GetCellCellAdhesionEnergy,
            " ", py::arg("pCellA"), py::arg("pCellB"))
        .def("GetCellBoundaryAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule_3::*)(::CellPtr)) &DifferentialAdhesionPottsUpdateRule_3::GetCellBoundaryAdhesionEnergy,
            " ", py::arg("pCell"))
        .def("GetLabelledCellLabelledCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule_3::*)()) &DifferentialAdhesionPottsUpdateRule_3::GetLabelledCellLabelledCellAdhesionEnergyParameter,
            " ")
        .def("GetLabelledCellCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule_3::*)()) &DifferentialAdhesionPottsUpdateRule_3::GetLabelledCellCellAdhesionEnergyParameter,
            " ")
        .def("GetLabelledCellBoundaryAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule_3::*)()) &DifferentialAdhesionPottsUpdateRule_3::GetLabelledCellBoundaryAdhesionEnergyParameter,
            " ")
        .def("SetLabelledCellLabelledCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule_3::*)(double)) &DifferentialAdhesionPottsUpdateRule_3::SetLabelledCellLabelledCellAdhesionEnergyParameter,
            " ", py::arg("labelledCellLabelledCellAdhesionEnergyParameter"))
        .def("SetLabelledCellCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule_3::*)(double)) &DifferentialAdhesionPottsUpdateRule_3::SetLabelledCellCellAdhesionEnergyParameter,
            " ", py::arg("labelledCellCellAdhesionEnergyParameter"))
        .def("SetLabelledCellBoundaryAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule_3::*)(double)) &DifferentialAdhesionPottsUpdateRule_3::SetLabelledCellBoundaryAdhesionEnergyParameter,
            " ", py::arg("labelledCellBoundaryAdhesionEnergyParameter"))
        .def("OutputUpdateRuleParameters",
            (void(DifferentialAdhesionPottsUpdateRule_3::*)(::out_stream &)) &DifferentialAdhesionPottsUpdateRule_3::OutputUpdateRuleParameters,
            " ", py::arg("rParamsFile"))
    ;
}
