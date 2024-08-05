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
#include "Toroidal2dMesh.hpp"

#include "Toroidal2dMesh.cppwg.hpp"

namespace py = pybind11;
typedef Toroidal2dMesh Toroidal2dMesh;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;

class Toroidal2dMesh_Overrides : public Toroidal2dMesh{
    public:
    using Toroidal2dMesh::Toroidal2dMesh;
    void ReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            ReMesh,
                    rMap);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorFromAtoB(::boost::numeric::ublas::c_vector<double, 2> const & rLocation1, ::boost::numeric::ublas::c_vector<double, 2> const & rLocation2) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            Toroidal2dMesh,
            GetVectorFromAtoB,
                    rLocation1,
        rLocation2);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> point, bool concreteMove) override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            SetNode,
                    index,
        point,
        concreteMove);
    }
    double GetWidth(unsigned int const & rDimension) const  override {
        PYBIND11_OVERRIDE(
            double,
            Toroidal2dMesh,
            GetWidth,
                    rDimension);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            Toroidal2dMesh,
            AddNode,
                    pNewNode);
    }
    void RefreshMesh() override {
        PYBIND11_OVERRIDE(
            void,
            Toroidal2dMesh,
            RefreshMesh,
            );
    }

};
void register_Toroidal2dMesh_class(py::module &m){
py::class_<Toroidal2dMesh , Toroidal2dMesh_Overrides , boost::shared_ptr<Toroidal2dMesh >  , MutableMesh<2, 2>  >(m, "Toroidal2dMesh")
        .def(py::init<double, double >(), py::arg("width"), py::arg("depth"))
        .def(py::init<double, double, ::std::vector<Node<2> *> >(), py::arg("width"), py::arg("depth"), py::arg("nodes"))
        .def(
            "ReMesh",
            (void(Toroidal2dMesh::*)(::NodeMap &)) &Toroidal2dMesh::ReMesh,
            " " , py::arg("rMap") )
        .def(
            "GetVectorFromAtoB",
            (::boost::numeric::ublas::c_vector<double, 2>(Toroidal2dMesh::*)(::boost::numeric::ublas::c_vector<double, 2> const &, ::boost::numeric::ublas::c_vector<double, 2> const &)) &Toroidal2dMesh::GetVectorFromAtoB,
            " " , py::arg("rLocation1"), py::arg("rLocation2") )
        .def(
            "SetNode",
            (void(Toroidal2dMesh::*)(unsigned int, ::ChastePoint<2>, bool)) &Toroidal2dMesh::SetNode,
            " " , py::arg("index"), py::arg("point"), py::arg("concreteMove") )
        .def(
            "GetWidth",
            (double(Toroidal2dMesh::*)(unsigned int const &) const ) &Toroidal2dMesh::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "AddNode",
            (unsigned int(Toroidal2dMesh::*)(::Node<2> *)) &Toroidal2dMesh::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "RefreshMesh",
            (void(Toroidal2dMesh::*)()) &Toroidal2dMesh::RefreshMesh,
            " "  )
    ;
}
