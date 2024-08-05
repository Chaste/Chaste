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
#include "NodeAttributes.hpp"

#include "NodeAttributes2.cppwg.hpp"

namespace py = pybind11;
typedef NodeAttributes<2 > NodeAttributes2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_NodeAttributes2_class(py::module &m){
py::class_<NodeAttributes2  , boost::shared_ptr<NodeAttributes2 >   >(m, "NodeAttributes2")
        .def(py::init< >())
        .def(
            "rGetAttributes",
            (::std::vector<double> &(NodeAttributes2::*)()) &NodeAttributes2::rGetAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "AddAttribute",
            (void(NodeAttributes2::*)(double)) &NodeAttributes2::AddAttribute,
            " " , py::arg("attribute") )
        .def(
            "GetRegion",
            (unsigned int(NodeAttributes2::*)()) &NodeAttributes2::GetRegion,
            " "  )
        .def(
            "SetRegion",
            (void(NodeAttributes2::*)(unsigned int)) &NodeAttributes2::SetRegion,
            " " , py::arg("region") )
        .def(
            "rGetAppliedForce",
            (::boost::numeric::ublas::c_vector<double, 2> &(NodeAttributes2::*)()) &NodeAttributes2::rGetAppliedForce,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "AddAppliedForceContribution",
            (void(NodeAttributes2::*)(::boost::numeric::ublas::c_vector<double, 2> const &)) &NodeAttributes2::AddAppliedForceContribution,
            " " , py::arg("rForceContribution") )
        .def(
            "ClearAppliedForce",
            (void(NodeAttributes2::*)()) &NodeAttributes2::ClearAppliedForce,
            " "  )
        .def(
            "AddNeighbour",
            (void(NodeAttributes2::*)(unsigned int)) &NodeAttributes2::AddNeighbour,
            " " , py::arg("index") )
        .def(
            "ClearNeighbours",
            (void(NodeAttributes2::*)()) &NodeAttributes2::ClearNeighbours,
            " "  )
        .def(
            "RemoveDuplicateNeighbours",
            (void(NodeAttributes2::*)()) &NodeAttributes2::RemoveDuplicateNeighbours,
            " "  )
        .def(
            "NeighboursIsEmpty",
            (bool(NodeAttributes2::*)()) &NodeAttributes2::NeighboursIsEmpty,
            " "  )
        .def(
            "SetNeighboursSetUp",
            (void(NodeAttributes2::*)(bool)) &NodeAttributes2::SetNeighboursSetUp,
            " " , py::arg("flag") )
        .def(
            "GetNeighboursSetUp",
            (bool(NodeAttributes2::*)()) &NodeAttributes2::GetNeighboursSetUp,
            " "  )
        .def(
            "rGetNeighbours",
            (::std::vector<unsigned int> &(NodeAttributes2::*)()) &NodeAttributes2::rGetNeighbours,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "IsParticle",
            (bool(NodeAttributes2::*)()) &NodeAttributes2::IsParticle,
            " "  )
        .def(
            "SetIsParticle",
            (void(NodeAttributes2::*)(bool)) &NodeAttributes2::SetIsParticle,
            " " , py::arg("isParticle") )
        .def(
            "GetRadius",
            (double(NodeAttributes2::*)()) &NodeAttributes2::GetRadius,
            " "  )
        .def(
            "SetRadius",
            (void(NodeAttributes2::*)(double)) &NodeAttributes2::SetRadius,
            " " , py::arg("radius") )
    ;
}
