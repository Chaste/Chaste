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
#include "Edge.hpp"

#include "Edge2.cppwg.hpp"

namespace py = pybind11;
typedef Edge<2 > Edge2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_Edge2_class(py::module &m){
py::class_<Edge2  , boost::shared_ptr<Edge2 >   >(m, "Edge2")
        .def(py::init<unsigned int >(), py::arg("index"))
        .def(py::init<unsigned int, ::Node<2> *, ::Node<2> * >(), py::arg("index"), py::arg("pNodeA"), py::arg("pNodeB"))
        .def_static(
            "GenerateMapIndex",
            (::std::pair<unsigned int, unsigned int>(*)(unsigned int, unsigned int)) &Edge2::GenerateMapIndex,
            " " , py::arg("index1"), py::arg("index2") )
        .def(
            "MarkAsDeleted",
            (void(Edge2::*)()) &Edge2::MarkAsDeleted,
            " "  )
        .def(
            "IsDeleted",
            (bool(Edge2::*)()) &Edge2::IsDeleted,
            " "  )
        .def(
            "SetIndex",
            (void(Edge2::*)(unsigned int)) &Edge2::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetIndex",
            (unsigned int(Edge2::*)() const ) &Edge2::GetIndex,
            " "  )
        .def(
            "GetMapIndex",
            (::std::pair<unsigned int, unsigned int>(Edge2::*)()) &Edge2::GetMapIndex,
            " "  )
        .def(
            "RemoveNodes",
            (void(Edge2::*)()) &Edge2::RemoveNodes,
            " "  )
        .def(
            "SetNodes",
            (void(Edge2::*)(::Node<2> *, ::Node<2> *)) &Edge2::SetNodes,
            " " , py::arg("pNodeA"), py::arg("pNodeB") )
        .def(
            "ReplaceNode",
            (void(Edge2::*)(::Node<2> *, ::Node<2> *)) &Edge2::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "GetNode",
            (::Node<2> *(Edge2::*)(unsigned int) const ) &Edge2::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(Edge2::*)()) &Edge2::GetNumNodes,
            " "  )
        .def(
            "ContainsNode",
            (bool(Edge2::*)(::Node<2> *) const ) &Edge2::ContainsNode,
            " " , py::arg("pNode") )
        .def(
            "rGetCentreLocation",
            (::boost::numeric::ublas::c_vector<double, 2>(Edge2::*)()) &Edge2::rGetCentreLocation,
            " "  )
        .def(
            "rGetLength",
            (double(Edge2::*)()) &Edge2::rGetLength,
            " "  )
        .def(
            "GetOtherElements",
            (::std::set<unsigned int>(Edge2::*)(unsigned int)) &Edge2::GetOtherElements,
            " " , py::arg("elementIndex") )
        .def(
            "AddElement",
            (void(Edge2::*)(unsigned int)) &Edge2::AddElement,
            " " , py::arg("elementIndex") )
        .def(
            "RemoveElement",
            (void(Edge2::*)(unsigned int)) &Edge2::RemoveElement,
            " " , py::arg("elementIndex") )
        .def(
            "GetNeighbouringElementIndices",
            (::std::set<unsigned int>(Edge2::*)()) &Edge2::GetNeighbouringElementIndices,
            " "  )
        .def(
            "GetNumElements",
            (unsigned int(Edge2::*)()) &Edge2::GetNumElements,
            " "  )
        .def(
            "IsBoundaryEdge",
            (bool(Edge2::*)() const ) &Edge2::IsBoundaryEdge,
            " "  )
    ;
}
