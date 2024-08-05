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

#include "ImmersedBoundarySimulationModifier2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundarySimulationModifier<2 > ImmersedBoundarySimulationModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundarySimulationModifier2_Overrides : public ImmersedBoundarySimulationModifier2{
    public:
    using ImmersedBoundarySimulationModifier2::ImmersedBoundarySimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySimulationModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundarySimulationModifier2_class(py::module &m){
py::class_<ImmersedBoundarySimulationModifier2 , ImmersedBoundarySimulationModifier2_Overrides , boost::shared_ptr<ImmersedBoundarySimulationModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "ImmersedBoundarySimulationModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ImmersedBoundarySimulationModifier2::*)(::AbstractCellPopulation<2> &)) &ImmersedBoundarySimulationModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ImmersedBoundarySimulationModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ImmersedBoundarySimulationModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ImmersedBoundarySimulationModifier2::*)(::out_stream &)) &ImmersedBoundarySimulationModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetNodeNeighbourUpdateFrequency",
            (void(ImmersedBoundarySimulationModifier2::*)(unsigned int)) &ImmersedBoundarySimulationModifier2::SetNodeNeighbourUpdateFrequency,
            " " , py::arg("newFrequency") )
        .def(
            "GetNodeNeighbourUpdateFrequency",
            (unsigned int(ImmersedBoundarySimulationModifier2::*)()) &ImmersedBoundarySimulationModifier2::GetNodeNeighbourUpdateFrequency,
            " "  )
        .def(
            "AddImmersedBoundaryForce",
            (void(ImmersedBoundarySimulationModifier2::*)(::boost::shared_ptr<AbstractImmersedBoundaryForce<2>>)) &ImmersedBoundarySimulationModifier2::AddImmersedBoundaryForce,
            " " , py::arg("pForce") )
        .def(
            "AddNormalNoise",
            (void(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::AddNormalNoise,
            " "  )
        .def(
            "GetZeroFieldSums",
            (bool(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetZeroFieldSums,
            " "  )
        .def(
            "SetZeroFieldSums",
            (void(ImmersedBoundarySimulationModifier2::*)(bool)) &ImmersedBoundarySimulationModifier2::SetZeroFieldSums,
            " " , py::arg("zeroFieldSums") )
        .def(
            "SetReynoldsNumber",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetReynoldsNumber,
            " " , py::arg("reynoldsNumber") )
        .def(
            "GetReynoldsNumber",
            (double(ImmersedBoundarySimulationModifier2::*)()) &ImmersedBoundarySimulationModifier2::GetReynoldsNumber,
            " "  )
        .def(
            "GetAdditiveNormalNoise",
            (bool(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetAdditiveNormalNoise,
            " "  )
        .def(
            "SetAdditiveNormalNoise",
            (void(ImmersedBoundarySimulationModifier2::*)(bool)) &ImmersedBoundarySimulationModifier2::SetAdditiveNormalNoise,
            " " , py::arg("additiveNormalNoise") )
        .def(
            "GetNoiseStrength",
            (double(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseStrength,
            " "  )
        .def(
            "SetNoiseStrength",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetNoiseStrength,
            " " , py::arg("noiseStrength") )
        .def(
            "GetNoiseSkip",
            (unsigned int(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseSkip,
            " "  )
        .def(
            "SetNoiseSkip",
            (void(ImmersedBoundarySimulationModifier2::*)(unsigned int)) &ImmersedBoundarySimulationModifier2::SetNoiseSkip,
            " " , py::arg("noiseSkip") )
        .def(
            "GetNoiseLengthScale",
            (double(ImmersedBoundarySimulationModifier2::*)() const ) &ImmersedBoundarySimulationModifier2::GetNoiseLengthScale,
            " "  )
        .def(
            "SetNoiseLengthScale",
            (void(ImmersedBoundarySimulationModifier2::*)(double)) &ImmersedBoundarySimulationModifier2::SetNoiseLengthScale,
            " " , py::arg("noiseLengthScale") )
    ;
}
