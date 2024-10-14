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
#include "PythonUblasObjectConverters.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"

#include "TargetAreaLinearGrowthModifier_3.cppwg.hpp"

namespace py = pybind11;
typedef TargetAreaLinearGrowthModifier<3> TargetAreaLinearGrowthModifier_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TargetAreaLinearGrowthModifier_3_Overrides : public TargetAreaLinearGrowthModifier_3
{
public:
    using TargetAreaLinearGrowthModifier_3::TargetAreaLinearGrowthModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier_3,
            UpdateTargetAreaOfCell,
            pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier_3,
            OutputSimulationModifierParameters,
            rParamsFile);
    }
};

void register_TargetAreaLinearGrowthModifier_3_class(py::module &m)
{
    py::class_<TargetAreaLinearGrowthModifier_3, TargetAreaLinearGrowthModifier_3_Overrides, boost::shared_ptr<TargetAreaLinearGrowthModifier_3>, AbstractTargetAreaModifier<3>>(m, "TargetAreaLinearGrowthModifier_3")
        .def(py::init<>())
        .def("UpdateTargetAreaOfCell",
            (void(TargetAreaLinearGrowthModifier_3::*)(::CellPtr const)) &TargetAreaLinearGrowthModifier_3::UpdateTargetAreaOfCell,
            " ", py::arg("pCell"))
        .def("GetAgeToStartGrowing",
            (double(TargetAreaLinearGrowthModifier_3::*)()) &TargetAreaLinearGrowthModifier_3::GetAgeToStartGrowing,
            " ")
        .def("SetAgeToStartGrowing",
            (void(TargetAreaLinearGrowthModifier_3::*)(double)) &TargetAreaLinearGrowthModifier_3::SetAgeToStartGrowing,
            " ", py::arg("ageToStartGrowing"))
        .def("GetGrowthRate",
            (double(TargetAreaLinearGrowthModifier_3::*)()) &TargetAreaLinearGrowthModifier_3::GetGrowthRate,
            " ")
        .def("SetGrowthRate",
            (void(TargetAreaLinearGrowthModifier_3::*)(double)) &TargetAreaLinearGrowthModifier_3::SetGrowthRate,
            " ", py::arg("growthRate"))
        .def("OutputSimulationModifierParameters",
            (void(TargetAreaLinearGrowthModifier_3::*)(::out_stream &)) &TargetAreaLinearGrowthModifier_3::OutputSimulationModifierParameters,
            " ", py::arg("rParamsFile"))
    ;
}
