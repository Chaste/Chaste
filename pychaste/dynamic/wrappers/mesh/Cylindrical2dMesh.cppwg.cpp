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
#include "Cylindrical2dMesh.hpp"

#include "Cylindrical2dMesh.cppwg.hpp"

namespace py = pybind11;
typedef Cylindrical2dMesh Cylindrical2dMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class Cylindrical2dMesh_Overrides : public Cylindrical2dMesh{
    public:
    using Cylindrical2dMesh::Cylindrical2dMesh;
    void ReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dMesh,
            ReMesh,
                    rMap);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Cylindrical2dMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            Cylindrical2dMesh,
            SetNode,
                    index,
        point,
        concreteMove);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Cylindrical2dMesh,
            GetWidth,
                    rDimension);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            Cylindrical2dMesh,
            AddNode,
                    pNewNode);
    }

};
void register_Cylindrical2dMesh_class(py::module &m){
py::class_<Cylindrical2dMesh , Cylindrical2dMesh_Overrides , boost::shared_ptr<Cylindrical2dMesh >  , MutableMesh<2, 2>  >(m, "Cylindrical2dMesh")
        .def(py::init<double >(), py::arg("width"))
        .def(py::init<double, ::std::vector<Node<2> *> >(), py::arg("width"), py::arg("nodes"))
        .def(
            "ReMesh",
            (void(Cylindrical2dMesh::*)(::NodeMap &)) &Cylindrical2dMesh::ReMesh,
            " " , py::arg("rMap") )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Cylindrical2dMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Cylindrical2dMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(Cylindrical2dMesh::*)(unsigned int, ::ChastePoint<2>, bool)) &Cylindrical2dMesh::SetNode,
            " " , py::arg("index"), py::arg("point"), py::arg("concreteMove") )
        .def(
            "GetWidth",
            (double(Cylindrical2dMesh::*)(unsigned int const &) const ) &Cylindrical2dMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "SetHaloScalingFactor",
            (void(Cylindrical2dMesh::*)(double)) &Cylindrical2dMesh::SetHaloScalingFactor,
            " " , py::arg("haloScalingFactor") )
        .def(
            "GetHaloScalingFactor",
            (double(Cylindrical2dMesh::*)() const ) &Cylindrical2dMesh::GetHaloScalingFactor,
            " "  )
        .def(
            "SetHaloOffset",
            (void(Cylindrical2dMesh::*)(double)) &Cylindrical2dMesh::SetHaloOffset,
            " " , py::arg("haloOffset") )
        .def(
            "GetHaloOffset",
            (double(Cylindrical2dMesh::*)() const ) &Cylindrical2dMesh::GetHaloOffset,
            " "  )
        .def(
            "AddNode",
            (unsigned int(Cylindrical2dMesh::*)(::Node<2> *)) &Cylindrical2dMesh::AddNode,
            " " , py::arg("pNewNode") )
    ;
}
