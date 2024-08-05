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
#include "CellwiseSourceEllipticPde.hpp"

#include "CellwiseSourceEllipticPde3.cppwg.hpp"

namespace py = pybind11;
typedef CellwiseSourceEllipticPde<3 > CellwiseSourceEllipticPde3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class CellwiseSourceEllipticPde3_Overrides : public CellwiseSourceEllipticPde3{
    public:
    using CellwiseSourceEllipticPde3::CellwiseSourceEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde3,
            ComputeConstantInUSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde3,
            ComputeLinearInUCoeffInSourceTerm,
                    rX,
        pElement);
    }
    double ComputeLinearInUCoeffInSourceTermAtNode(::Node<3> const & rNode) override {
        PYBIND11_OVERRIDE(
            double,
            CellwiseSourceEllipticPde3,
            ComputeLinearInUCoeffInSourceTermAtNode,
                    rNode);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            CellwiseSourceEllipticPde3,
            ComputeDiffusionTerm,
                    rX);
    }

};
void register_CellwiseSourceEllipticPde3_class(py::module &m){
py::class_<CellwiseSourceEllipticPde3 , CellwiseSourceEllipticPde3_Overrides , boost::shared_ptr<CellwiseSourceEllipticPde3 >  , AbstractLinearEllipticPde<3, 3>  >(m, "CellwiseSourceEllipticPde3")
        .def(py::init<::AbstractCellPopulation<3> &, double >(), py::arg("rCellPopulation"), py::arg("sourceCoefficient") = 0.)
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<3> const &(CellwiseSourceEllipticPde3::*)() const ) &CellwiseSourceEllipticPde3::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetCoefficient",
            (double(CellwiseSourceEllipticPde3::*)() const ) &CellwiseSourceEllipticPde3::GetCoefficient,
            " "  )
        .def(
            "ComputeConstantInUSourceTerm",
            (double(CellwiseSourceEllipticPde3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &CellwiseSourceEllipticPde3::ComputeConstantInUSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTerm",
            (double(CellwiseSourceEllipticPde3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &CellwiseSourceEllipticPde3::ComputeLinearInUCoeffInSourceTerm,
            " " , py::arg("rX"), py::arg("pElement") )
        .def(
            "ComputeLinearInUCoeffInSourceTermAtNode",
            (double(CellwiseSourceEllipticPde3::*)(::Node<3> const &)) &CellwiseSourceEllipticPde3::ComputeLinearInUCoeffInSourceTermAtNode,
            " " , py::arg("rNode") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(CellwiseSourceEllipticPde3::*)(::ChastePoint<3> const &)) &CellwiseSourceEllipticPde3::ComputeDiffusionTerm,
            " " , py::arg("rX") )
    ;
}
