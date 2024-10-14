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
#include "PythonUblasObjectConverters.hpp"
#include "PeriodicNodesOnlyMesh.hpp"

#include "PeriodicNodesOnlyMesh_2.cppwg.hpp"

namespace py = pybind11;
typedef PeriodicNodesOnlyMesh<2> PeriodicNodesOnlyMesh_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class PeriodicNodesOnlyMesh_2_Overrides : public PeriodicNodesOnlyMesh_2
{
public:
    using PeriodicNodesOnlyMesh_2::PeriodicNodesOnlyMesh;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            PeriodicNodesOnlyMesh_2,
            GetVectorFromAtoB,
            rLocation1,
            rLocation2);
    }
    double GetWidth(unsigned int const & rDimension) const override
    {
        PYBIND11_OVERRIDE(
            double,
            PeriodicNodesOnlyMesh_2,
            GetWidth,
            rDimension);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> point, bool concreteMove) override
    {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh_2,
            SetNode,
            nodeIndex,
            point,
            concreteMove);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            PeriodicNodesOnlyMesh_2,
            AddNode,
            pNewNode);
    }
    void RefreshMesh() override
    {
        PYBIND11_OVERRIDE(
            void,
            PeriodicNodesOnlyMesh_2,
            RefreshMesh,
            );
    }
};

void register_PeriodicNodesOnlyMesh_2_class(py::module &m)
{
    py::class_<PeriodicNodesOnlyMesh_2, PeriodicNodesOnlyMesh_2_Overrides, boost::shared_ptr<PeriodicNodesOnlyMesh_2>, NodesOnlyMesh<2>>(m, "PeriodicNodesOnlyMesh_2")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2>>(), py::arg("width"))
        .def("GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(PeriodicNodesOnlyMesh_2::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &PeriodicNodesOnlyMesh_2::GetVectorFromAtoB,
            " ", py::arg("rLocation1"), py::arg("rLocation2"))
        .def("GetWidth",
            (double(PeriodicNodesOnlyMesh_2::*)(unsigned int const &) const) &PeriodicNodesOnlyMesh_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("GetPeriodicWidths",
            (::boost::numeric::ublas::c_vector<double, 2>(PeriodicNodesOnlyMesh_2::*)() const) &PeriodicNodesOnlyMesh_2::GetPeriodicWidths,
            " ")
        .def("SetNode",
            (void(PeriodicNodesOnlyMesh_2::*)(unsigned int, ::ChastePoint<2>, bool)) &PeriodicNodesOnlyMesh_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("point"), py::arg("concreteMove") = false)
        .def("AddNode",
            (unsigned int(PeriodicNodesOnlyMesh_2::*)(::Node<2> *)) &PeriodicNodesOnlyMesh_2::AddNode,
            " ", py::arg("pNewNode"))
        .def("RefreshMesh",
            (void(PeriodicNodesOnlyMesh_2::*)()) &PeriodicNodesOnlyMesh_2::RefreshMesh,
            " ")
    ;
}
