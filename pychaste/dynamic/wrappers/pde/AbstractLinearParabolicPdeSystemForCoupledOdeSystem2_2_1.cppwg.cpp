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
#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearParabolicPdeSystemForCoupledOdeSystem<2,2,1 > AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_Overrides : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1{
    public:
    using AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::AbstractLinearParabolicPdeSystemForCoupledOdeSystem;
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeDuDtCoefficientFunction,
                    rX,
        pdeIndex);
    }
    double ComputeSourceTerm(::ChastePoint<2> const & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeSourceTerm,
                    rX,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    double ComputeSourceTermAtNode(::Node<2> const & rNode, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::std::vector<double> & rOdeSolution, unsigned int pdeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeSourceTermAtNode,
                    rNode,
        rU,
        rOdeSolution,
        pdeIndex);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, unsigned int pdeIndex, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1,
            ComputeDiffusionTerm,
                    rX,
        pdeIndex,
        pElement);
    }

};
void register_AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_class(py::module &m){
py::class_<AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1 , AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1_Overrides , boost::shared_ptr<AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1 >   >(m, "AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1")
        .def(py::init< >())
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTerm",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeSourceTermAtNode",
            (double(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::Node<2> const &, ::boost::numeric::ublas::c_vector<double, 1> &, ::std::vector<double> &, unsigned int)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeSourceTermAtNode,
            " " , py::arg("rNode"), py::arg("rU"), py::arg("rOdeSolution"), py::arg("pdeIndex") )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::*)(::ChastePoint<2> const &, unsigned int, ::Element<2, 2> *)) &AbstractLinearParabolicPdeSystemForCoupledOdeSystem2_2_1::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pdeIndex"), py::arg("pElement") = __null )
    ;
}
