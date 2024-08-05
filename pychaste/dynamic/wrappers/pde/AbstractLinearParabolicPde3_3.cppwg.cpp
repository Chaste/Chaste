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
#include "AbstractLinearParabolicPde.hpp"

#include "AbstractLinearParabolicPde3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPde<3,3 > AbstractLinearParabolicPde3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AbstractLinearParabolicPde3_3_Overrides : public AbstractLinearParabolicPde3_3{
    public:
    using AbstractLinearParabolicPde3_3::AbstractLinearParabolicPde;
    double ComputeDuDtCoefficientFunction(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeDuDtCoefficientFunction,
                    rX);
    }
    double ComputeSourceTerm(::ChastePoint<3> const & rX, double u, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    double ComputeSourceTermAtNode(::Node<3> const & rNode, double u) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPde3_3,
            ComputeSourceTermAtNode,
                    rNode,
        u);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AbstractLinearParabolicPde3_3,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }

};
void register_AbstractLinearParabolicPde3_3_class(py::module &m){
py::class_<AbstractLinearParabolicPde3_3 , AbstractLinearParabolicPde3_3_Overrides , boost::shared_ptr<AbstractLinearParabolicPde3_3 >  , AbstractLinearPde<3>  >(m, "AbstractLinearParabolicPde3_3")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &)) &AbstractLinearParabolicPde3_3::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &, double, ::Element<3, 3> *)) &AbstractLinearParabolicPde3_3::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = nullptr )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPde3_3::*)(::Node<3> const &, double)) &AbstractLinearParabolicPde3_3::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("u") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AbstractLinearParabolicPde3_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AbstractLinearParabolicPde3_3::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
    ;
}
