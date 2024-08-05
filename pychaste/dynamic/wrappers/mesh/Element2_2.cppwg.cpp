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
#include "Element.hpp"

#include "Element2_2.cppwg.hpp"

namespace py = pybind11;
typedef Element<2,2 > Element2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class Element2_2_Overrides : public Element2_2{
    public:
    using Element2_2::Element;
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            MarkAsDeleted,
            );
    }
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            UpdateNode,
                    rIndex,
        pNode);
    }

};
void register_Element2_2_class(py::module &m){
py::class_<Element2_2 , Element2_2_Overrides , boost::shared_ptr<Element2_2 >   >(m, "Element2_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const &, bool >(), py::arg("index"), py::arg("rNodes"), py::arg("registerWithNodes") = true)
        .def(py::init<::Element<2, 2> const &, unsigned int const >(), py::arg("rElement"), py::arg("index"))
        .def(
            "RegisterWithNodes",
            (void(Element2_2::*)()) &Element2_2::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(Element2_2::*)()) &Element2_2::MarkAsDeleted,
            " "  )
        .def(
            "UpdateNode",
            (void(Element2_2::*)(unsigned int const &, ::Node<2> *)) &Element2_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ResetIndex",
            (void(Element2_2::*)(unsigned int)) &Element2_2::ResetIndex,
            " " , py::arg("index") )
        .def(
            "CalculateCircumsphere",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::boost::numeric::ublas::c_matrix<double, 2, 2> &, ::boost::numeric::ublas::c_matrix<double, 2, 2> &)) &Element2_2::CalculateCircumsphere,
            " " , py::arg("rJacobian"), py::arg("rInverseJacobian") )
        .def(
            "CalculateQuality",
            (double(Element2_2::*)()) &Element2_2::CalculateQuality,
            " "  )
        .def(
            "CalculateMinMaxEdgeLengths",
            (::boost::numeric::ublas::c_vector<double, 2>(Element2_2::*)()) &Element2_2::CalculateMinMaxEdgeLengths,
            " "  )
        .def(
            "CalculateInterpolationWeights",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateInterpolationWeights,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateInterpolationWeightsWithProjection",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateInterpolationWeightsWithProjection,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateXi",
            (::boost::numeric::ublas::c_vector<double, 2>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateXi,
            " " , py::arg("rTestPoint") )
        .def(
            "IncludesPoint",
            (bool(Element2_2::*)(::ChastePoint<2> const &, bool)) &Element2_2::IncludesPoint,
            " " , py::arg("rTestPoint"), py::arg("strict") = false )
    ;
}
