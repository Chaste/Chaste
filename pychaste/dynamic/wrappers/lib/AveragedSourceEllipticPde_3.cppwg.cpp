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
#include "AveragedSourceEllipticPde.hpp"

#include "AveragedSourceEllipticPde_3.cppwg.hpp"

namespace py = pybind11;
typedef AveragedSourceEllipticPde<3> AveragedSourceEllipticPde_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AveragedSourceEllipticPde_3_Overrides : public AveragedSourceEllipticPde_3
{
public:
    using AveragedSourceEllipticPde_3::AveragedSourceEllipticPde;
    void SetupSourceTerms(::TetrahedralMesh<3, 3> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            AveragedSourceEllipticPde_3,
            SetupSourceTerms,
            rCoarseMesh,
            pCellPdeElementMap);
    }
    double ComputeConstantInUSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceEllipticPde_3,
            ComputeConstantInUSourceTerm,
            rX,
            pElement);
    }
    double ComputeLinearInUCoeffInSourceTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceEllipticPde_3,
            ComputeLinearInUCoeffInSourceTerm,
            rX,
            pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AveragedSourceEllipticPde_3,
            ComputeDiffusionTerm,
            rX);
    }
};

void register_AveragedSourceEllipticPde_3_class(py::module &m)
{
    py::class_<AveragedSourceEllipticPde_3, AveragedSourceEllipticPde_3_Overrides, boost::shared_ptr<AveragedSourceEllipticPde_3>, AbstractLinearEllipticPde<3, 3>>(m, "AveragedSourceEllipticPde_3")
        .def(py::init<::AbstractCellPopulation<3> &, double, double>(), py::arg("rCellPopulation"), py::arg("sourceCoefficient") = 0., py::arg("diffusionCoefficient") = 1.)
        .def("rGetCellPopulation",
            (::AbstractCellPopulation<3> const &(AveragedSourceEllipticPde_3::*)() const) &AveragedSourceEllipticPde_3::rGetCellPopulation,
            " ", py::return_value_policy::reference_internal)
        .def("GetCoefficient",
            (double(AveragedSourceEllipticPde_3::*)() const) &AveragedSourceEllipticPde_3::GetCoefficient,
            " ")
        .def("SetupSourceTerms",
            (void(AveragedSourceEllipticPde_3::*)(::TetrahedralMesh<3, 3> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AveragedSourceEllipticPde_3::SetupSourceTerms,
            " ", py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr)
        .def("ComputeConstantInUSourceTerm",
            (double(AveragedSourceEllipticPde_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AveragedSourceEllipticPde_3::ComputeConstantInUSourceTerm,
            " ", py::arg("rX"), py::arg("pElement"))
        .def("ComputeLinearInUCoeffInSourceTerm",
            (double(AveragedSourceEllipticPde_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AveragedSourceEllipticPde_3::ComputeLinearInUCoeffInSourceTerm,
            " ", py::arg("rX"), py::arg("pElement"))
        .def("ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AveragedSourceEllipticPde_3::*)(::ChastePoint<3> const &)) &AveragedSourceEllipticPde_3::ComputeDiffusionTerm,
            " ", py::arg("rX"))
        .def("GetUptakeRateForElement",
            (double(AveragedSourceEllipticPde_3::*)(unsigned int)) &AveragedSourceEllipticPde_3::GetUptakeRateForElement,
            " ", py::arg("elementIndex"))
    ;
}
