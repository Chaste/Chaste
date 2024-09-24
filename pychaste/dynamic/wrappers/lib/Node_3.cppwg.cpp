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

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Node.hpp"

#include "Node_3.cppwg.hpp"

namespace py = pybind11;
typedef Node<3> Node_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_Node_3_class(py::module &m)
{
    py::class_<Node_3, boost::shared_ptr<Node_3>>(m, "Node_3")
        .def(py::init<unsigned int, ::ChastePoint<3>, bool>(), py::arg("index"), py::arg("point"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, ::std::vector<double>, bool>(), py::arg("index"), py::arg("coords"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, ::boost::numeric::ublas::c_vector<double, 3>, bool>(), py::arg("index"), py::arg("location"), py::arg("isBoundaryNode") = false)
        .def(py::init<unsigned int, bool, double, double, double>(), py::arg("index"), py::arg("isBoundaryNode") = false, py::arg("v1") = 0, py::arg("v2") = 0, py::arg("v3") = 0)
        .def(py::init<unsigned int, double *, bool>(), py::arg("index"), py::arg("location"), py::arg("isBoundaryNode") = false)
        .def("SetPoint",
            (void(Node_3::*)(::ChastePoint<3>)) &Node_3::SetPoint,
            " ", py::arg("point"))
        .def("SetIndex",
            (void(Node_3::*)(unsigned int)) &Node_3::SetIndex,
            " ", py::arg("index"))
        .def("AddNodeAttribute",
            (void(Node_3::*)(double)) &Node_3::AddNodeAttribute,
            " ", py::arg("attribute"))
        .def("SetAsBoundaryNode",
            (void(Node_3::*)(bool)) &Node_3::SetAsBoundaryNode,
            " ", py::arg("value") = true)
        .def("GetPoint",
            (::ChastePoint<3>(Node_3::*)() const) &Node_3::GetPoint,
            " ")
        .def("rGetLocation",
            (::boost::numeric::ublas::c_vector<double, 3> const &(Node_3::*)() const) &Node_3::rGetLocation,
            " ", py::return_value_policy::reference_internal)
        .def("rGetModifiableLocation",
            (::boost::numeric::ublas::c_vector<double, 3> &(Node_3::*)()) &Node_3::rGetModifiableLocation,
            " ", py::return_value_policy::reference_internal)
        .def("GetIndex",
            (unsigned int(Node_3::*)() const) &Node_3::GetIndex,
            " ")
        .def("IsBoundaryNode",
            (bool(Node_3::*)() const) &Node_3::IsBoundaryNode,
            " ")
        .def("AddElement",
            (void(Node_3::*)(unsigned int)) &Node_3::AddElement,
            " ", py::arg("index"))
        .def("RemoveElement",
            (void(Node_3::*)(unsigned int)) &Node_3::RemoveElement,
            " ", py::arg("index"))
        .def("RemoveBoundaryElement",
            (void(Node_3::*)(unsigned int)) &Node_3::RemoveBoundaryElement,
            " ", py::arg("index"))
        .def("AddBoundaryElement",
            (void(Node_3::*)(unsigned int)) &Node_3::AddBoundaryElement,
            " ", py::arg("index"))
        .def("AddNeighbour",
            (void(Node_3::*)(unsigned int)) &Node_3::AddNeighbour,
            " ", py::arg("index"))
        .def("ClearNeighbours",
            (void(Node_3::*)()) &Node_3::ClearNeighbours,
            " ")
        .def("RemoveDuplicateNeighbours",
            (void(Node_3::*)()) &Node_3::RemoveDuplicateNeighbours,
            " ")
        .def("NeighboursIsEmpty",
            (bool(Node_3::*)()) &Node_3::NeighboursIsEmpty,
            " ")
        .def("SetNeighboursSetUp",
            (void(Node_3::*)(bool)) &Node_3::SetNeighboursSetUp,
            " ", py::arg("flag"))
        .def("GetNeighboursSetUp",
            (bool(Node_3::*)()) &Node_3::GetNeighboursSetUp,
            " ")
        .def("rGetNeighbours",
            (::std::vector<unsigned int> &(Node_3::*)()) &Node_3::rGetNeighbours,
            " ", py::return_value_policy::reference_internal)
        .def("rGetContainingElementIndices",
            (::std::set<unsigned int> &(Node_3::*)()) &Node_3::rGetContainingElementIndices,
            " ", py::return_value_policy::reference_internal)
        .def("rGetNodeAttributes",
            (::std::vector<double> &(Node_3::*)()) &Node_3::rGetNodeAttributes,
            " ", py::return_value_policy::reference_internal)
        .def("GetNumNodeAttributes",
            (unsigned int(Node_3::*)()) &Node_3::GetNumNodeAttributes,
            " ")
        .def("HasNodeAttributes",
            (bool(Node_3::*)()) &Node_3::HasNodeAttributes,
            " ")
        .def("rGetContainingBoundaryElementIndices",
            (::std::set<unsigned int> &(Node_3::*)()) &Node_3::rGetContainingBoundaryElementIndices,
            " ", py::return_value_policy::reference_internal)
        .def("GetNumContainingElements",
            (unsigned int(Node_3::*)() const) &Node_3::GetNumContainingElements,
            " ")
        .def("GetNumBoundaryElements",
            (unsigned int(Node_3::*)() const) &Node_3::GetNumBoundaryElements,
            " ")
        .def("rGetAppliedForce",
            (::boost::numeric::ublas::c_vector<double, 3> &(Node_3::*)()) &Node_3::rGetAppliedForce,
            " ", py::return_value_policy::reference_internal)
        .def("ClearAppliedForce",
            (void(Node_3::*)()) &Node_3::ClearAppliedForce,
            " ")
        .def("AddAppliedForceContribution",
            (void(Node_3::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &Node_3::AddAppliedForceContribution,
            " ", py::arg("rForceContribution"))
        .def("IsParticle",
            (bool(Node_3::*)()) &Node_3::IsParticle,
            " ")
        .def("SetIsParticle",
            (void(Node_3::*)(bool)) &Node_3::SetIsParticle,
            " ", py::arg("isParticle"))
        .def("GetRadius",
            (double(Node_3::*)()) &Node_3::GetRadius,
            " ")
        .def("SetRadius",
            (void(Node_3::*)(double)) &Node_3::SetRadius,
            " ", py::arg("radius"))
        .def("MarkAsDeleted",
            (void(Node_3::*)()) &Node_3::MarkAsDeleted,
            " ")
        .def("IsDeleted",
            (bool(Node_3::*)() const) &Node_3::IsDeleted,
            " ")
        .def("MarkAsInternal",
            (void(Node_3::*)()) &Node_3::MarkAsInternal,
            " ")
        .def("IsInternal",
            (bool(Node_3::*)() const) &Node_3::IsInternal,
            " ")
        .def("SetRegion",
            (void(Node_3::*)(unsigned int)) &Node_3::SetRegion,
            " ", py::arg("region"))
        .def("GetRegion",
            (unsigned int(Node_3::*)() const) &Node_3::GetRegion,
            " ")
        .def("ContainingElementsBegin",
            (::Node<3>::ContainingElementIterator(Node_3::*)() const) &Node_3::ContainingElementsBegin,
            " ")
        .def("ContainingElementsEnd",
            (::Node<3>::ContainingElementIterator(Node_3::*)() const) &Node_3::ContainingElementsEnd,
            " ")
        .def("ContainingBoundaryElementsBegin",
            (::Node<3>::ContainingBoundaryElementIterator(Node_3::*)() const) &Node_3::ContainingBoundaryElementsBegin,
            " ")
        .def("ContainingBoundaryElementsEnd",
            (::Node<3>::ContainingBoundaryElementIterator(Node_3::*)() const) &Node_3::ContainingBoundaryElementsEnd,
            " ")
    ;
}
