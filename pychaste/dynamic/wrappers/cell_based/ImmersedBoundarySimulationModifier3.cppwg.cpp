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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"

#include "ImmersedBoundarySimulationModifier3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundarySimulationModifier<3 > ImmersedBoundarySimulationModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundarySimulationModifier3_Overrides : public ImmersedBoundarySimulationModifier3{
    public:
    using ImmersedBoundarySimulationModifier3::ImmersedBoundarySimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundarySimulationModifier3_class(py::module &m){
py::class_<ImmersedBoundarySimulationModifier3 , ImmersedBoundarySimulationModifier3_Overrides , boost::shared_ptr<ImmersedBoundarySimulationModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "ImmersedBoundarySimulationModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ImmersedBoundarySimulationModifier3::*)(::AbstractCellPopulation<3> &)) &ImmersedBoundarySimulationModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ImmersedBoundarySimulationModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &ImmersedBoundarySimulationModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ImmersedBoundarySimulationModifier3::*)(::out_stream &)) &ImmersedBoundarySimulationModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetNodeNeighbourUpdateFrequency",
            (void(ImmersedBoundarySimulationModifier3::*)(unsigned int)) &ImmersedBoundarySimulationModifier3::SetNodeNeighbourUpdateFrequency,
            " " , py::arg("newFrequency") )
        .def(
            "GetNodeNeighbourUpdateFrequency",
            (unsigned int(ImmersedBoundarySimulationModifier3::*)()) &ImmersedBoundarySimulationModifier3::GetNodeNeighbourUpdateFrequency,
            " "  )
        .def(
            "AddImmersedBoundaryForce",
            (void(ImmersedBoundarySimulationModifier3::*)(::boost::shared_ptr<AbstractImmersedBoundaryForce<3>>)) &ImmersedBoundarySimulationModifier3::AddImmersedBoundaryForce,
            " " , py::arg("pForce") )
        .def(
            "AddNormalNoise",
            (void(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::AddNormalNoise,
            " "  )
        .def(
            "GetZeroFieldSums",
            (bool(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::GetZeroFieldSums,
            " "  )
        .def(
            "SetZeroFieldSums",
            (void(ImmersedBoundarySimulationModifier3::*)(bool)) &ImmersedBoundarySimulationModifier3::SetZeroFieldSums,
            " " , py::arg("zeroFieldSums") )
        .def(
            "SetReynoldsNumber",
            (void(ImmersedBoundarySimulationModifier3::*)(double)) &ImmersedBoundarySimulationModifier3::SetReynoldsNumber,
            " " , py::arg("reynoldsNumber") )
        .def(
            "GetReynoldsNumber",
            (double(ImmersedBoundarySimulationModifier3::*)()) &ImmersedBoundarySimulationModifier3::GetReynoldsNumber,
            " "  )
        .def(
            "GetAdditiveNormalNoise",
            (bool(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::GetAdditiveNormalNoise,
            " "  )
        .def(
            "SetAdditiveNormalNoise",
            (void(ImmersedBoundarySimulationModifier3::*)(bool)) &ImmersedBoundarySimulationModifier3::SetAdditiveNormalNoise,
            " " , py::arg("additiveNormalNoise") )
        .def(
            "GetNoiseStrength",
            (double(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::GetNoiseStrength,
            " "  )
        .def(
            "SetNoiseStrength",
            (void(ImmersedBoundarySimulationModifier3::*)(double)) &ImmersedBoundarySimulationModifier3::SetNoiseStrength,
            " " , py::arg("noiseStrength") )
        .def(
            "GetNoiseSkip",
            (unsigned int(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::GetNoiseSkip,
            " "  )
        .def(
            "SetNoiseSkip",
            (void(ImmersedBoundarySimulationModifier3::*)(unsigned int)) &ImmersedBoundarySimulationModifier3::SetNoiseSkip,
            " " , py::arg("noiseSkip") )
        .def(
            "GetNoiseLengthScale",
            (double(ImmersedBoundarySimulationModifier3::*)() const ) &ImmersedBoundarySimulationModifier3::GetNoiseLengthScale,
            " "  )
        .def(
            "SetNoiseLengthScale",
            (void(ImmersedBoundarySimulationModifier3::*)(double)) &ImmersedBoundarySimulationModifier3::SetNoiseLengthScale,
            " " , py::arg("noiseLengthScale") )
    ;
}
