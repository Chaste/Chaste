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
#include "AbstractLinearEllipticPde.hpp"

#include "AbstractLinearEllipticPde_2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearEllipticPde<2, 2> AbstractLinearEllipticPde_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 2, 2> _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_;

class AbstractLinearEllipticPde_2_2_Overrides : public AbstractLinearEllipticPde_2_2
{
public:
    using AbstractLinearEllipticPde_2_2::AbstractLinearEllipticPde;
    double ComputeConstantInUSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde_2_2,
            ComputeConstantInUSourceTerm,
            rX,
            pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<2> const & rX, ::Element<2, 2> * pElement) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractLinearEllipticPde_2_2,
            ComputeLinearInUCoeffInSourceTerm,
            rX,
            pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 2, 2> ComputeDiffusionTerm(::ChastePoint<2> const & rX) override
    {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_matrix_lt_double_2_2_gt_,
            AbstractLinearEllipticPde_2_2,
            ComputeDiffusionTerm,
            rX);
    }
    double ComputeConstantInUSourceTermAtNode(::Node<2> const & rNode) override
    {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde_2_2,
            ComputeConstantInUSourceTermAtNode,
            rNode);
    }
    double ComputeLinearInUCoeffInSourceTermAtNode(::Node<2> const & rNode) override
    {
        PYBIND11_OVERRIDE(
            double,
            AbstractLinearEllipticPde_2_2,
            ComputeLinearInUCoeffInSourceTermAtNode,
            rNode);
    }
};

void register_AbstractLinearEllipticPde_2_2_class(py::module &m)
{
    py::class_<AbstractLinearEllipticPde_2_2, AbstractLinearEllipticPde_2_2_Overrides, boost::shared_ptr<AbstractLinearEllipticPde_2_2>, AbstractLinearPde<2>>(m, "AbstractLinearEllipticPde_2_2")
        .def(py::init<>())
        .def("ComputeConstantInUSourceTerm",
            (double(AbstractLinearEllipticPde_2_2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AbstractLinearEllipticPde_2_2::ComputeConstantInUSourceTerm,
            " ", py::arg("rX"), py::arg("pElement"))
        .def("ComputeLinearInUCoeffInSourceTerm",
            (double(AbstractLinearEllipticPde_2_2::*)(::ChastePoint<2> const &, ::Element<2, 2> *)) &AbstractLinearEllipticPde_2_2::ComputeLinearInUCoeffInSourceTerm,
            " ", py::arg("rX"), py::arg("pElement"))
        .def("ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 2, 2>(AbstractLinearEllipticPde_2_2::*)(::ChastePoint<2> const &)) &AbstractLinearEllipticPde_2_2::ComputeDiffusionTerm,
            " ", py::arg("rX"))
        .def("ComputeConstantInUSourceTermAtNode",
            (double(AbstractLinearEllipticPde_2_2::*)(::Node<2> const &)) &AbstractLinearEllipticPde_2_2::ComputeConstantInUSourceTermAtNode,
            " ", py::arg("rNode"))
        .def("ComputeLinearInUCoeffInSourceTermAtNode",
            (double(AbstractLinearEllipticPde_2_2::*)(::Node<2> const &)) &AbstractLinearEllipticPde_2_2::ComputeLinearInUCoeffInSourceTermAtNode,
            " ", py::arg("rNode"))
    ;
}
