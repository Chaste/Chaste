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
#include "UniformSourceParabolicPde.hpp"

#include "UniformSourceParabolicPde2.cppwg.hpp"

namespace py = pybind11;
typedef UniformSourceParabolicPde<2 > UniformSourceParabolicPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class UniformSourceParabolicPde2_Overrides : public UniformSourceParabolicPde2{
    public:
    using UniformSourceParabolicPde2::UniformSourceParabolicPde;
    double ComputeSourceTerm(::ChastePoint<2> const & rX, double u, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde2,
            ComputeSourceTerm,
                    rX,
        u,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            UniformSourceParabolicPde2,
            ComputeDiffusionTerm,
                    rX,
        pElement);
    }
    double ComputeDuDtCoefficientFunction(::ChastePoint<2> const & rX) override {
        PYBIND11_OVERRIDE(
            double,
            UniformSourceParabolicPde2,
            ComputeDuDtCoefficientFunction,
                    rX);
    }

};
void register_UniformSourceParabolicPde2_class(py::module &m){
py::class_<UniformSourceParabolicPde2 , UniformSourceParabolicPde2_Overrides , boost::shared_ptr<UniformSourceParabolicPde2 >  , AbstractLinearParabolicPde<2>  >(m, "UniformSourceParabolicPde2")
        .def(py::init<double >(), py::arg("sourceCoefficient") = 0.)
        .def(
            "GetCoefficient",
            (double(UniformSourceParabolicPde2::*)() const ) &UniformSourceParabolicPde2::GetCoefficient,
            " "  )
        .def(
            "ComputeSourceTerm",
            (double(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &, double, ::Element<2, 2> *)) &UniformSourceParabolicPde2::ComputeSourceTerm,
            " " , py::arg("rX"), py::arg("u"), py::arg("pElement") = __null )
        .def(
            "ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &UniformSourceParabolicPde2::ComputeDiffusionTerm,
            " " , py::arg("rX"), py::arg("pElement") = __null )
        .def(
            "ComputeDuDtCoefficientFunction",
            (double(UniformSourceParabolicPde2::*)(::ChastePoint<2> const &)) &UniformSourceParabolicPde2::ComputeDuDtCoefficientFunction,
            " " , py::arg("rX") )
    ;
}
