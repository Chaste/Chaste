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
#include "ImmersedBoundarySvgWriter.hpp"

#include "ImmersedBoundarySvgWriter2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundarySvgWriter<2 > ImmersedBoundarySvgWriter2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundarySvgWriter2_Overrides : public ImmersedBoundarySvgWriter2{
    public:
    using ImmersedBoundarySvgWriter2::ImmersedBoundarySvgWriter;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySvgWriter2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySvgWriter2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundarySvgWriter2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundarySvgWriter2_class(py::module &m){
py::class_<ImmersedBoundarySvgWriter2 , ImmersedBoundarySvgWriter2_Overrides , boost::shared_ptr<ImmersedBoundarySvgWriter2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "ImmersedBoundarySvgWriter2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ImmersedBoundarySvgWriter2::*)(::AbstractCellPopulation<2> &)) &ImmersedBoundarySvgWriter2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ImmersedBoundarySvgWriter2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ImmersedBoundarySvgWriter2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "AddPointToSvgFile",
            (void(ImmersedBoundarySvgWriter2::*)(::out_stream &, ::boost::numeric::ublas::c_vector<double, 2>, unsigned int, double)) &ImmersedBoundarySvgWriter2::AddPointToSvgFile,
            " " , py::arg("rSvgFile"), py::arg("location"), py::arg("region"), py::arg("rad") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ImmersedBoundarySvgWriter2::*)(::out_stream &)) &ImmersedBoundarySvgWriter2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetSamplingMultiple",
            (unsigned int(ImmersedBoundarySvgWriter2::*)() const ) &ImmersedBoundarySvgWriter2::GetSamplingMultiple,
            " "  )
        .def(
            "SetSamplingMultiple",
            (void(ImmersedBoundarySvgWriter2::*)(unsigned int)) &ImmersedBoundarySvgWriter2::SetSamplingMultiple,
            " " , py::arg("samplingMultiple") )
        .def(
            "GetSvgSize",
            (double(ImmersedBoundarySvgWriter2::*)() const ) &ImmersedBoundarySvgWriter2::GetSvgSize,
            " "  )
        .def(
            "SetSvgSize",
            (void(ImmersedBoundarySvgWriter2::*)(double)) &ImmersedBoundarySvgWriter2::SetSvgSize,
            " " , py::arg("svgSize") )
    ;
}
