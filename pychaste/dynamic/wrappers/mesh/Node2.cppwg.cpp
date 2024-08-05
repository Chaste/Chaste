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
#include "Node.hpp"

#include "Node2.cppwg.hpp"

namespace py = pybind11;
typedef Node<2 > Node2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_Node2_class(py::module &m){
py::class_<Node2  , boost::shared_ptr<Node2 >   >(m, "Node2")
        .def(py::init<unsigned int, ::ChastePoint<2>, bool >(), py::arg("index"), py::arg("point"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, ::std::vector<double>, bool >(), py::arg("index"), py::arg("coords"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, ::boost::numeric::ublas::c_vector<double, 2>, bool >(), py::arg("index"), py::arg("location"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, bool, double, double, double >(), py::arg("index"), py::arg("isBoundaryNode") = false, py::arg("v1") = 0, py::arg("v2") = 0, py::arg("v3") = 0)
        .def(py::init<unsigned int, double *, bool >(), py::arg("index"), py::arg("location"), py::arg("isBoundaryNode") = false)
        .def(
            "SetPoint",
            (void(Node2::*)(::ChastePoint<2>)) &Node2::SetPoint,
            " " , py::arg("point") )
        .def(
            "SetIndex",
            (void(Node2::*)(unsigned int)) &Node2::SetIndex,
            " " , py::arg("index") )
        .def(
            "AddNodeAttribute",
            (void(Node2::*)(double)) &Node2::AddNodeAttribute,
            " " , py::arg("attribute") )
        .def(
            "SetAsBoundaryNode",
            (void(Node2::*)(bool)) &Node2::SetAsBoundaryNode,
            " " , py::arg("value") = true )
        .def(
            "GetPoint",
            (::ChastePoint<2>(Node2::*)() const ) &Node2::GetPoint,
            " "  )
        .def(
            "rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 2> const &(Node2::*)() const ) &Node2::rGetLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetModifiableLocation",
            (::boost::numeric::ublas::c_vector<double, 2> &(Node2::*)()) &Node2::rGetModifiableLocation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetIndex",
            (unsigned int(Node2::*)() const ) &Node2::GetIndex,
            " "  )
        .def(
            "IsBoundaryNode",
            (bool(Node2::*)() const ) &Node2::IsBoundaryNode,
            " "  )
        .def(
            "AddElement",
            (void(Node2::*)(unsigned int)) &Node2::AddElement,
            " " , py::arg("index") )
        .def(
            "RemoveElement",
            (void(Node2::*)(unsigned int)) &Node2::RemoveElement,
            " " , py::arg("index") )
        .def(
            "RemoveBoundaryElement",
            (void(Node2::*)(unsigned int)) &Node2::RemoveBoundaryElement,
            " " , py::arg("index") )
        .def(
            "AddBoundaryElement",
            (void(Node2::*)(unsigned int)) &Node2::AddBoundaryElement,
            " " , py::arg("index") )
        .def(
            "AddNeighbour",
            (void(Node2::*)(unsigned int)) &Node2::AddNeighbour,
            " " , py::arg("index") )
        .def(
            "ClearNeighbours",
            (void(Node2::*)()) &Node2::ClearNeighbours,
            " "  )
        .def(
            "RemoveDuplicateNeighbours",
            (void(Node2::*)()) &Node2::RemoveDuplicateNeighbours,
            " "  )
        .def(
            "NeighboursIsEmpty",
            (bool(Node2::*)()) &Node2::NeighboursIsEmpty,
            " "  )
        .def(
            "SetNeighboursSetUp",
            (void(Node2::*)(bool)) &Node2::SetNeighboursSetUp,
            " " , py::arg("flag") )
        .def(
            "GetNeighboursSetUp",
            (bool(Node2::*)()) &Node2::GetNeighboursSetUp,
            " "  )
        .def(
            "rGetNeighbours",
            (::std::vector<unsigned int> &(Node2::*)()) &Node2::rGetNeighbours,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetContainingElementIndices",
            (::std::set<unsigned int> &(Node2::*)()) &Node2::rGetContainingElementIndices,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetNodeAttributes",
            (::std::vector<double> &(Node2::*)()) &Node2::rGetNodeAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumNodeAttributes",
            (unsigned int(Node2::*)()) &Node2::GetNumNodeAttributes,
            " "  )
        .def(
            "HasNodeAttributes",
            (bool(Node2::*)()) &Node2::HasNodeAttributes,
            " "  )
        .def(
            "rGetContainingBoundaryElementIndices",
            (::std::set<unsigned int> &(Node2::*)()) &Node2::rGetContainingBoundaryElementIndices,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumContainingElements",
            (unsigned int(Node2::*)() const ) &Node2::GetNumContainingElements,
            " "  )
        .def(
            "GetNumBoundaryElements",
            (unsigned int(Node2::*)() const ) &Node2::GetNumBoundaryElements,
            " "  )
        .def(
            "rGetAppliedForce",
            (::boost::numeric::ublas::c_vector<double, 2> &(Node2::*)()) &Node2::rGetAppliedForce,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "ClearAppliedForce",
            (void(Node2::*)()) &Node2::ClearAppliedForce,
            " "  )
        .def(
            "AddAppliedForceContribution",
            (void(Node2::*)(::boost::numeric::ublas::c_vector<double, 2> const &)) &Node2::AddAppliedForceContribution,
            " " , py::arg("rForceContribution") )
        .def(
            "IsParticle",
            (bool(Node2::*)()) &Node2::IsParticle,
            " "  )
        .def(
            "SetIsParticle",
            (void(Node2::*)(bool)) &Node2::SetIsParticle,
            " " , py::arg("isParticle") )
        .def(
            "GetRadius",
            (double(Node2::*)()) &Node2::GetRadius,
            " "  )
        .def(
            "SetRadius",
            (void(Node2::*)(double)) &Node2::SetRadius,
            " " , py::arg("radius") )
        .def(
            "MarkAsDeleted",
            (void(Node2::*)()) &Node2::MarkAsDeleted,
            " "  )
        .def(
            "IsDeleted",
            (bool(Node2::*)() const ) &Node2::IsDeleted,
            " "  )
        .def(
            "MarkAsInternal",
            (void(Node2::*)()) &Node2::MarkAsInternal,
            " "  )
        .def(
            "IsInternal",
            (bool(Node2::*)() const ) &Node2::IsInternal,
            " "  )
        .def(
            "SetRegion",
            (void(Node2::*)(unsigned int)) &Node2::SetRegion,
            " " , py::arg("region") )
        .def(
            "GetRegion",
            (unsigned int(Node2::*)() const ) &Node2::GetRegion,
            " "  )
        .def(
            "ContainingElementsBegin",
            (::Node<2>::ContainingElementIterator(Node2::*)() const ) &Node2::ContainingElementsBegin,
            " "  )
        .def(
            "ContainingElementsEnd",
            (::Node<2>::ContainingElementIterator(Node2::*)() const ) &Node2::ContainingElementsEnd,
            " "  )
        .def(
            "ContainingBoundaryElementsBegin",
            (::Node<2>::ContainingBoundaryElementIterator(Node2::*)() const ) &Node2::ContainingBoundaryElementsBegin,
            " "  )
        .def(
            "ContainingBoundaryElementsEnd",
            (::Node<2>::ContainingBoundaryElementIterator(Node2::*)() const ) &Node2::ContainingBoundaryElementsEnd,
            " "  )
    ;
}
