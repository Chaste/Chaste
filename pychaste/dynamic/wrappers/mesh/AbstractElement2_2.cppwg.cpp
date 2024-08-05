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
#include "AbstractElement.hpp"

#include "AbstractElement2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractElement<2,2 > AbstractElement2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractElement2_2_Overrides : public AbstractElement2_2{
    public:
    using AbstractElement2_2::AbstractElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            UpdateNode,
                    rIndex,
        pNode);
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            MarkAsDeleted,
            );
    }
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            RegisterWithNodes,
            );
    }

};
void register_AbstractElement2_2_class(py::module &m){
py::class_<AbstractElement2_2 , AbstractElement2_2_Overrides , boost::shared_ptr<AbstractElement2_2 >   >(m, "AbstractElement2_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(py::init<unsigned int >(), py::arg("index") = ::INDEX_IS_NOT_USED)
        .def(
            "UpdateNode",
            (void(AbstractElement2_2::*)(unsigned int const &, ::Node<2> *)) &AbstractElement2_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ReplaceNode",
            (void(AbstractElement2_2::*)(::Node<2> *, ::Node<2> *)) &AbstractElement2_2::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "MarkAsDeleted",
            (void(AbstractElement2_2::*)()) &AbstractElement2_2::MarkAsDeleted,
            " "  )
        .def(
            "RegisterWithNodes",
            (void(AbstractElement2_2::*)()) &AbstractElement2_2::RegisterWithNodes,
            " "  )
        .def(
            "GetNodeLocation",
            (double(AbstractElement2_2::*)(unsigned int, unsigned int) const ) &AbstractElement2_2::GetNodeLocation,
            " " , py::arg("localIndex"), py::arg("dimension") )
        .def(
            "GetNodeLocation",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNodeLocation,
            " " , py::arg("localIndex") )
        .def(
            "GetNodeGlobalIndex",
            (unsigned int(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNodeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNode",
            (::Node<2> *(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNode,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetNumNodes,
            " "  )
        .def(
            "AddNode",
            (void(AbstractElement2_2::*)(::Node<2> *)) &AbstractElement2_2::AddNode,
            " " , py::arg("pNode") )
        .def(
            "IsDeleted",
            (bool(AbstractElement2_2::*)() const ) &AbstractElement2_2::IsDeleted,
            " "  )
        .def(
            "GetIndex",
            (unsigned int(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(AbstractElement2_2::*)(unsigned int)) &AbstractElement2_2::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetOwnership",
            (bool(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetOwnership,
            " "  )
        .def(
            "SetOwnership",
            (void(AbstractElement2_2::*)(bool)) &AbstractElement2_2::SetOwnership,
            " " , py::arg("ownership") )
        .def(
            "SetAttribute",
            (void(AbstractElement2_2::*)(double)) &AbstractElement2_2::SetAttribute,
            " " , py::arg("attribute") )
        .def(
            "GetAttribute",
            (double(AbstractElement2_2::*)()) &AbstractElement2_2::GetAttribute,
            " "  )
        .def(
            "GetUnsignedAttribute",
            (unsigned int(AbstractElement2_2::*)()) &AbstractElement2_2::GetUnsignedAttribute,
            " "  )
        .def(
            "AddElementAttribute",
            (void(AbstractElement2_2::*)(double)) &AbstractElement2_2::AddElementAttribute,
            " " , py::arg("attribute") )
        .def(
            "rGetElementAttributes",
            (::std::vector<double> &(AbstractElement2_2::*)()) &AbstractElement2_2::rGetElementAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumElementAttributes",
            (unsigned int(AbstractElement2_2::*)()) &AbstractElement2_2::GetNumElementAttributes,
            " "  )
    ;
}
