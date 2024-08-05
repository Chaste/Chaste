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
#include "ImmersedBoundaryMorseInteractionForce.hpp"

#include "ImmersedBoundaryMorseInteractionForce3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseInteractionForce<3 > ImmersedBoundaryMorseInteractionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseInteractionForce3_Overrides : public ImmersedBoundaryMorseInteractionForce3{
    public:
    using ImmersedBoundaryMorseInteractionForce3::ImmersedBoundaryMorseInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<3> *, Node<3> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce3,
            AddImmersedBoundaryForceContribution,
                    rNodePairs,
        rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce3,
            OutputImmersedBoundaryForceParameters,
                    rParamsFile);
    }

};
void register_ImmersedBoundaryMorseInteractionForce3_class(py::module &m){
py::class_<ImmersedBoundaryMorseInteractionForce3 , ImmersedBoundaryMorseInteractionForce3_Overrides , boost::shared_ptr<ImmersedBoundaryMorseInteractionForce3 >  , AbstractImmersedBoundaryForce<3>  >(m, "ImmersedBoundaryMorseInteractionForce3")
        .def(py::init< >())
        .def(
            "AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(::std::vector<std::pair<Node<3> *, Node<3> *>> &, ::ImmersedBoundaryCellPopulation<3> &)) &ImmersedBoundaryMorseInteractionForce3::AddImmersedBoundaryForceContribution,
            " " , py::arg("rNodePairs"), py::arg("rCellPopulation") )
        .def(
            "OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(::out_stream &)) &ImmersedBoundaryMorseInteractionForce3::OutputImmersedBoundaryForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetWellDepth",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetWellDepth,
            " "  )
        .def(
            "SetWellDepth",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetWellDepth,
            " " , py::arg("wellDepth") )
        .def(
            "GetRestLength",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetRestLength,
            " "  )
        .def(
            "SetRestLength",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetRestLength,
            " " , py::arg("restLength") )
        .def(
            "GetLaminaWellDepthMult",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetLaminaWellDepthMult,
            " "  )
        .def(
            "SetLaminaWellDepthMult",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetLaminaWellDepthMult,
            " " , py::arg("laminaWellDepthMult") )
        .def(
            "GetLaminaRestLengthMult",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetLaminaRestLengthMult,
            " "  )
        .def(
            "SetLaminaRestLengthMult",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetLaminaRestLengthMult,
            " " , py::arg("laminaRestLengthMult") )
        .def(
            "GetWellWidth",
            (double(ImmersedBoundaryMorseInteractionForce3::*)() const ) &ImmersedBoundaryMorseInteractionForce3::GetWellWidth,
            " "  )
        .def(
            "SetWellWidth",
            (void(ImmersedBoundaryMorseInteractionForce3::*)(double)) &ImmersedBoundaryMorseInteractionForce3::SetWellWidth,
            " " , py::arg("wellWidth") )
    ;
}
