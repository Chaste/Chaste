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
#include "ImmersedBoundaryMorseInteractionForce.hpp"

#include "ImmersedBoundaryMorseInteractionForce_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryMorseInteractionForce<2> ImmersedBoundaryMorseInteractionForce_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryMorseInteractionForce_2_Overrides : public ImmersedBoundaryMorseInteractionForce_2
{
public:
    using ImmersedBoundaryMorseInteractionForce_2::ImmersedBoundaryMorseInteractionForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce_2,
            AddImmersedBoundaryForceContribution,
            rNodePairs,
            rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryMorseInteractionForce_2,
            OutputImmersedBoundaryForceParameters,
            rParamsFile);
    }
};

void register_ImmersedBoundaryMorseInteractionForce_2_class(py::module &m)
{
    py::class_<ImmersedBoundaryMorseInteractionForce_2, ImmersedBoundaryMorseInteractionForce_2_Overrides, boost::shared_ptr<ImmersedBoundaryMorseInteractionForce_2>, AbstractImmersedBoundaryForce<2>>(m, "ImmersedBoundaryMorseInteractionForce_2")
        .def(py::init<>())
        .def("AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryMorseInteractionForce_2::AddImmersedBoundaryForceContribution,
            " ", py::arg("rNodePairs"), py::arg("rCellPopulation"))
        .def("OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(::out_stream &)) &ImmersedBoundaryMorseInteractionForce_2::OutputImmersedBoundaryForceParameters,
            " ", py::arg("rParamsFile"))
        .def("GetWellDepth",
            (double(ImmersedBoundaryMorseInteractionForce_2::*)() const) &ImmersedBoundaryMorseInteractionForce_2::GetWellDepth,
            " ")
        .def("SetWellDepth",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(double)) &ImmersedBoundaryMorseInteractionForce_2::SetWellDepth,
            " ", py::arg("wellDepth"))
        .def("GetRestLength",
            (double(ImmersedBoundaryMorseInteractionForce_2::*)() const) &ImmersedBoundaryMorseInteractionForce_2::GetRestLength,
            " ")
        .def("SetRestLength",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(double)) &ImmersedBoundaryMorseInteractionForce_2::SetRestLength,
            " ", py::arg("restLength"))
        .def("GetLaminaWellDepthMult",
            (double(ImmersedBoundaryMorseInteractionForce_2::*)() const) &ImmersedBoundaryMorseInteractionForce_2::GetLaminaWellDepthMult,
            " ")
        .def("SetLaminaWellDepthMult",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(double)) &ImmersedBoundaryMorseInteractionForce_2::SetLaminaWellDepthMult,
            " ", py::arg("laminaWellDepthMult"))
        .def("GetLaminaRestLengthMult",
            (double(ImmersedBoundaryMorseInteractionForce_2::*)() const) &ImmersedBoundaryMorseInteractionForce_2::GetLaminaRestLengthMult,
            " ")
        .def("SetLaminaRestLengthMult",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(double)) &ImmersedBoundaryMorseInteractionForce_2::SetLaminaRestLengthMult,
            " ", py::arg("laminaRestLengthMult"))
        .def("GetWellWidth",
            (double(ImmersedBoundaryMorseInteractionForce_2::*)() const) &ImmersedBoundaryMorseInteractionForce_2::GetWellWidth,
            " ")
        .def("SetWellWidth",
            (void(ImmersedBoundaryMorseInteractionForce_2::*)(double)) &ImmersedBoundaryMorseInteractionForce_2::SetWellWidth,
            " ", py::arg("wellWidth"))
    ;
}
