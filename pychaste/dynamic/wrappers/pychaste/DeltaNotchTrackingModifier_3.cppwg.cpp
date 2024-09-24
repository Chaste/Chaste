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

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchTrackingModifier.hpp"

#include "DeltaNotchTrackingModifier_3.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchTrackingModifier<3> DeltaNotchTrackingModifier_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchTrackingModifier_3_Overrides : public DeltaNotchTrackingModifier_3
{
public:
    using DeltaNotchTrackingModifier_3::DeltaNotchTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier_3,
            UpdateAtEndOfTimeStep,
            rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier_3,
            SetupSolve,
            rCellPopulation,
            outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier_3,
            OutputSimulationModifierParameters,
            rParamsFile);
    }
};

void register_DeltaNotchTrackingModifier_3_class(py::module &m)
{
    py::class_<DeltaNotchTrackingModifier_3, DeltaNotchTrackingModifier_3_Overrides, boost::shared_ptr<DeltaNotchTrackingModifier_3>, AbstractCellBasedSimulationModifier<3>>(m, "DeltaNotchTrackingModifier_3")
        .def(py::init<>())
        .def("UpdateAtEndOfTimeStep",
            (void(DeltaNotchTrackingModifier_3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchTrackingModifier_3::UpdateAtEndOfTimeStep,
            " ", py::arg("rCellPopulation"))
        .def("SetupSolve",
            (void(DeltaNotchTrackingModifier_3::*)(::AbstractCellPopulation<3> &, ::std::string)) &DeltaNotchTrackingModifier_3::SetupSolve,
            " ", py::arg("rCellPopulation"), py::arg("outputDirectory"))
        .def("UpdateCellData",
            (void(DeltaNotchTrackingModifier_3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchTrackingModifier_3::UpdateCellData,
            " ", py::arg("rCellPopulation"))
        .def("OutputSimulationModifierParameters",
            (void(DeltaNotchTrackingModifier_3::*)(::out_stream &)) &DeltaNotchTrackingModifier_3::OutputSimulationModifierParameters,
            " ", py::arg("rParamsFile"))
    ;
}
