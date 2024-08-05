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
#include "AbstractTargetAreaModifier.hpp"

#include "AbstractTargetAreaModifier3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTargetAreaModifier<3 > AbstractTargetAreaModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractTargetAreaModifier3_Overrides : public AbstractTargetAreaModifier3{
    public:
    using AbstractTargetAreaModifier3::AbstractTargetAreaModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractTargetAreaModifier3,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractTargetAreaModifier3_class(py::module &m){
py::class_<AbstractTargetAreaModifier3 , AbstractTargetAreaModifier3_Overrides , boost::shared_ptr<AbstractTargetAreaModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "AbstractTargetAreaModifier3")
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractTargetAreaModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &AbstractTargetAreaModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "GetReferenceTargetArea",
            (double(AbstractTargetAreaModifier3::*)()) &AbstractTargetAreaModifier3::GetReferenceTargetArea,
            " "  )
        .def(
            "SetReferenceTargetArea",
            (void(AbstractTargetAreaModifier3::*)(double)) &AbstractTargetAreaModifier3::SetReferenceTargetArea,
            " " , py::arg("referenceTargetArea") )
        .def(
            "UpdateTargetAreas",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractTargetAreaModifier3::UpdateTargetAreas,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateTargetAreaOfCell",
            (void(AbstractTargetAreaModifier3::*)(::CellPtr const)) &AbstractTargetAreaModifier3::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractTargetAreaModifier3::*)(::out_stream &)) &AbstractTargetAreaModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
