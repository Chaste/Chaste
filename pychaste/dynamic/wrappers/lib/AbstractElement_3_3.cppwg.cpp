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
#include "AbstractElement.hpp"

#include "AbstractElement_3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractElement<3, 3> AbstractElement_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractElement_3_3_Overrides : public AbstractElement_3_3
{
public:
    using AbstractElement_3_3::AbstractElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<3> * pNode) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement_3_3,
            UpdateNode,
            rIndex,
            pNode);
    }
    void MarkAsDeleted() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement_3_3,
            MarkAsDeleted,
            );
    }
    void RegisterWithNodes() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement_3_3,
            RegisterWithNodes,
            );
    }
};

void register_AbstractElement_3_3_class(py::module &m)
{
    py::class_<AbstractElement_3_3, AbstractElement_3_3_Overrides, boost::shared_ptr<AbstractElement_3_3>>(m, "AbstractElement_3_3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const &>(), py::arg("index"), py::arg("rNodes"))
        .def(py::init<unsigned int>(), py::arg("index") = ::INDEX_IS_NOT_USED)
        .def("UpdateNode",
            (void(AbstractElement_3_3::*)(unsigned int const &, ::Node<3> *)) &AbstractElement_3_3::UpdateNode,
            " ", py::arg("rIndex"), py::arg("pNode"))
        .def("ReplaceNode",
            (void(AbstractElement_3_3::*)(::Node<3> *, ::Node<3> *)) &AbstractElement_3_3::ReplaceNode,
            " ", py::arg("pOldNode"), py::arg("pNewNode"))
        .def("MarkAsDeleted",
            (void(AbstractElement_3_3::*)()) &AbstractElement_3_3::MarkAsDeleted,
            " ")
        .def("RegisterWithNodes",
            (void(AbstractElement_3_3::*)()) &AbstractElement_3_3::RegisterWithNodes,
            " ")
        .def("GetNodeLocation",
            (double(AbstractElement_3_3::*)(unsigned int, unsigned int) const) &AbstractElement_3_3::GetNodeLocation,
            " ", py::arg("localIndex"), py::arg("dimension"))
        .def("GetNodeLocation",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractElement_3_3::*)(unsigned int) const) &AbstractElement_3_3::GetNodeLocation,
            " ", py::arg("localIndex"))
        .def("GetNodeGlobalIndex",
            (unsigned int(AbstractElement_3_3::*)(unsigned int) const) &AbstractElement_3_3::GetNodeGlobalIndex,
            " ", py::arg("localIndex"))
        .def("GetNode",
            (::Node<3> *(AbstractElement_3_3::*)(unsigned int) const) &AbstractElement_3_3::GetNode,
            " ", py::arg("localIndex"), py::return_value_policy::reference)
        .def("GetNumNodes",
            (unsigned int(AbstractElement_3_3::*)() const) &AbstractElement_3_3::GetNumNodes,
            " ")
        .def("AddNode",
            (void(AbstractElement_3_3::*)(::Node<3> *)) &AbstractElement_3_3::AddNode,
            " ", py::arg("pNode"))
        .def("IsDeleted",
            (bool(AbstractElement_3_3::*)() const) &AbstractElement_3_3::IsDeleted,
            " ")
        .def("GetIndex",
            (unsigned int(AbstractElement_3_3::*)() const) &AbstractElement_3_3::GetIndex,
            " ")
        .def("SetIndex",
            (void(AbstractElement_3_3::*)(unsigned int)) &AbstractElement_3_3::SetIndex,
            " ", py::arg("index"))
        .def("GetOwnership",
            (bool(AbstractElement_3_3::*)() const) &AbstractElement_3_3::GetOwnership,
            " ")
        .def("SetOwnership",
            (void(AbstractElement_3_3::*)(bool)) &AbstractElement_3_3::SetOwnership,
            " ", py::arg("ownership"))
        .def("SetAttribute",
            (void(AbstractElement_3_3::*)(double)) &AbstractElement_3_3::SetAttribute,
            " ", py::arg("attribute"))
        .def("GetAttribute",
            (double(AbstractElement_3_3::*)()) &AbstractElement_3_3::GetAttribute,
            " ")
        .def("GetUnsignedAttribute",
            (unsigned int(AbstractElement_3_3::*)()) &AbstractElement_3_3::GetUnsignedAttribute,
            " ")
        .def("AddElementAttribute",
            (void(AbstractElement_3_3::*)(double)) &AbstractElement_3_3::AddElementAttribute,
            " ", py::arg("attribute"))
        .def("rGetElementAttributes",
            (::std::vector<double> &(AbstractElement_3_3::*)()) &AbstractElement_3_3::rGetElementAttributes,
            " ", py::return_value_policy::reference_internal)
        .def("GetNumElementAttributes",
            (unsigned int(AbstractElement_3_3::*)()) &AbstractElement_3_3::GetNumElementAttributes,
            " ")
    ;
}
