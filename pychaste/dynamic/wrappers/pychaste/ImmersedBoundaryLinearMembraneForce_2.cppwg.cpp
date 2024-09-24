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
#include "ImmersedBoundaryLinearMembraneForce.hpp"

#include "ImmersedBoundaryLinearMembraneForce_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryLinearMembraneForce<2> ImmersedBoundaryLinearMembraneForce_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryLinearMembraneForce_2_Overrides : public ImmersedBoundaryLinearMembraneForce_2
{
public:
    using ImmersedBoundaryLinearMembraneForce_2::ImmersedBoundaryLinearMembraneForce;
    void AddImmersedBoundaryForceContribution(::std::vector<std::pair<Node<2> *, Node<2> *>> & rNodePairs, ::ImmersedBoundaryCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce_2,
            AddImmersedBoundaryForceContribution,
            rNodePairs,
            rCellPopulation);
    }
    void OutputImmersedBoundaryForceParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryLinearMembraneForce_2,
            OutputImmersedBoundaryForceParameters,
            rParamsFile);
    }
};

void register_ImmersedBoundaryLinearMembraneForce_2_class(py::module &m)
{
    py::class_<ImmersedBoundaryLinearMembraneForce_2, ImmersedBoundaryLinearMembraneForce_2_Overrides, boost::shared_ptr<ImmersedBoundaryLinearMembraneForce_2>, AbstractImmersedBoundaryForce<2>>(m, "ImmersedBoundaryLinearMembraneForce_2")
        .def(py::init<>())
        .def("AddImmersedBoundaryForceContribution",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(::std::vector<std::pair<Node<2> *, Node<2> *>> &, ::ImmersedBoundaryCellPopulation<2> &)) &ImmersedBoundaryLinearMembraneForce_2::AddImmersedBoundaryForceContribution,
            " ", py::arg("rNodePairs"), py::arg("rCellPopulation"))
        .def("OutputImmersedBoundaryForceParameters",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(::out_stream &)) &ImmersedBoundaryLinearMembraneForce_2::OutputImmersedBoundaryForceParameters,
            " ", py::arg("rParamsFile"))
        .def("GetElementSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce_2::*)() const) &ImmersedBoundaryLinearMembraneForce_2::GetElementSpringConst,
            " ")
        .def("SetElementSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(double)) &ImmersedBoundaryLinearMembraneForce_2::SetElementSpringConst,
            " ", py::arg("elementSpringConst"))
        .def("GetElementRestLength",
            (double(ImmersedBoundaryLinearMembraneForce_2::*)() const) &ImmersedBoundaryLinearMembraneForce_2::GetElementRestLength,
            " ")
        .def("SetElementRestLength",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(double)) &ImmersedBoundaryLinearMembraneForce_2::SetElementRestLength,
            " ", py::arg("elementRestLength"))
        .def("GetLaminaSpringConst",
            (double(ImmersedBoundaryLinearMembraneForce_2::*)() const) &ImmersedBoundaryLinearMembraneForce_2::GetLaminaSpringConst,
            " ")
        .def("SetLaminaSpringConst",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(double)) &ImmersedBoundaryLinearMembraneForce_2::SetLaminaSpringConst,
            " ", py::arg("laminaSpringConst"))
        .def("GetLaminaRestLength",
            (double(ImmersedBoundaryLinearMembraneForce_2::*)() const) &ImmersedBoundaryLinearMembraneForce_2::GetLaminaRestLength,
            " ")
        .def("SetLaminaRestLength",
            (void(ImmersedBoundaryLinearMembraneForce_2::*)(double)) &ImmersedBoundaryLinearMembraneForce_2::SetLaminaRestLength,
            " ", py::arg("laminaRestLength"))
    ;
}
