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

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractNonlinearEllipticPde.hpp"

#include "AbstractNonlinearEllipticPde_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNonlinearEllipticPde<2> AbstractNonlinearEllipticPde_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractNonlinearEllipticPde_2_Overrides : public AbstractNonlinearEllipticPde_2
{
public:
    using AbstractNonlinearEllipticPde_2::AbstractNonlinearEllipticPde;
    double ComputeLinearSourceTerm(::ChastePoint<2> const & rX) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde_2,
            ComputeLinearSourceTerm,
            rX);
    }
    double ComputeNonlinearSourceTerm(::ChastePoint<2> const & rX, double u) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde_2,
            ComputeNonlinearSourceTerm,
            rX,
            u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX, double u) override
    {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractNonlinearEllipticPde_2,
            ComputeDiffusionTerm,
            rX,
            u);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTermPrime(::ChastePoint<2> const & rX, double u) override
    {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractNonlinearEllipticPde_2,
            ComputeDiffusionTermPrime,
            rX,
            u);
    }
    double ComputeNonlinearSourceTermPrime(::ChastePoint<2> const & rX, double u) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractNonlinearEllipticPde_2,
            ComputeNonlinearSourceTermPrime,
            rX,
            u);
    }
};

void register_AbstractNonlinearEllipticPde_2_class(py::module &m)
{
    py::class_<AbstractNonlinearEllipticPde_2, AbstractNonlinearEllipticPde_2_Overrides, boost::shared_ptr<AbstractNonlinearEllipticPde_2>>(m, "AbstractNonlinearEllipticPde_2")
        .def(py::init<>())
        .def("ComputeLinearSourceTerm",
            (double(AbstractNonlinearEllipticPde_2::*)(::ChastePoint<2> const &)) &AbstractNonlinearEllipticPde_2::ComputeLinearSourceTerm,
            " ", py::arg("rX"))
        .def("ComputeNonlinearSourceTerm",
            (double(AbstractNonlinearEllipticPde_2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde_2::ComputeNonlinearSourceTerm,
            " ", py::arg("rX"), py::arg("u"))
        .def("ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractNonlinearEllipticPde_2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde_2::ComputeDiffusionTerm,
            " ", py::arg("rX"), py::arg("u"))
        .def("ComputeDiffusionTermPrime",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractNonlinearEllipticPde_2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde_2::ComputeDiffusionTermPrime,
            " ", py::arg("rX"), py::arg("u"))
        .def("ComputeNonlinearSourceTermPrime",
            (double(AbstractNonlinearEllipticPde_2::*)(::ChastePoint<2> const &, double)) &AbstractNonlinearEllipticPde_2::ComputeNonlinearSourceTermPrime,
            " ", py::arg("rX"), py::arg("u"))
    ;
}
