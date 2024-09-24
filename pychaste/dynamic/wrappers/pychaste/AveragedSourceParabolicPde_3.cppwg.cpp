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
#include "AveragedSourceParabolicPde.hpp"

#include "AveragedSourceParabolicPde_3.cppwg.hpp"

namespace py = pybind11;
typedef AveragedSourceParabolicPde<3> AveragedSourceParabolicPde_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_matrix<double, 3, 3> _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_;

class AveragedSourceParabolicPde_3_Overrides : public AveragedSourceParabolicPde_3
{
public:
    using AveragedSourceParabolicPde_3::AveragedSourceParabolicPde;
    void SetupSourceTerms(::TetrahedralMesh<3, 3> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            AveragedSourceParabolicPde_3,
            SetupSourceTerms,
            rCoarseMesh,
            pCellPdeElementMap);
    }
    double ComputeDuDtCoefficientFunction(::ChastePoint<3> const & rX) override
    {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde_3,
            ComputeDuDtCoefficientFunction,
            rX);
    }
    double ComputeSourceTerm(::ChastePoint<3> const & rX, double u, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde_3,
            ComputeSourceTerm,
            rX,
            u,
            pElement);
    }
    double ComputeSourceTermAtNode(::Node<3> const & rNode, double u) override
    {
        PYBIND11_OVERRIDE(
            double,
            AveragedSourceParabolicPde_3,
            ComputeSourceTermAtNode,
            rNode,
            u);
    }
    ::boost::numeric::ublas::c_matrix<double, 3, 3> ComputeDiffusionTerm(::ChastePoint<3> const & rX, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_3_3_gt_,
            AveragedSourceParabolicPde_3,
            ComputeDiffusionTerm,
            rX,
            pElement);
    }
};

void register_AveragedSourceParabolicPde_3_class(py::module &m)
{
    py::class_<AveragedSourceParabolicPde_3, AveragedSourceParabolicPde_3_Overrides, boost::shared_ptr<AveragedSourceParabolicPde_3>, AbstractLinearParabolicPde<3>>(m, "AveragedSourceParabolicPde_3")
        .def(py::init<::AbstractCellPopulation<3> &, double, double, double>(), py::arg("rCellPopulation"), py::arg("duDtCoefficient") = 1., py::arg("diffusionCoefficient") = 1., py::arg("sourceCoefficient") = 0.)
        .def("rGetCellPopulation",
            (::AbstractCellPopulation<3> const &(AveragedSourceParabolicPde_3::*)() const) &AveragedSourceParabolicPde_3::rGetCellPopulation,
            " ", py::return_value_policy::reference_internal)
        .def("SetupSourceTerms",
            (void(AveragedSourceParabolicPde_3::*)(::TetrahedralMesh<3, 3> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &AveragedSourceParabolicPde_3::SetupSourceTerms,
            " ", py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr)
        .def("ComputeDuDtCoefficientFunction",
            (double(AveragedSourceParabolicPde_3::*)(::ChastePoint<3> const &)) &AveragedSourceParabolicPde_3::ComputeDuDtCoefficientFunction,
            " ", py::arg("rX"))
        .def("ComputeSourceTerm",
            (double(AveragedSourceParabolicPde_3::*)(::ChastePoint<3> const &, double, ::Element<3, 3> *)) &AveragedSourceParabolicPde_3::ComputeSourceTerm,
            " ", py::arg("rX"), py::arg("u"), py::arg("pElement") = __null)
        .def("ComputeSourceTermAtNode",
            (double(AveragedSourceParabolicPde_3::*)(::Node<3> const &, double)) &AveragedSourceParabolicPde_3::ComputeSourceTermAtNode,
            " ", py::arg("rNode"), py::arg("u"))
        .def("ComputeDiffusionTerm",
            (::boost::numeric::ublas::c_matrix<double, 3, 3>(AveragedSourceParabolicPde_3::*)(::ChastePoint<3> const &, ::Element<3, 3> *)) &AveragedSourceParabolicPde_3::ComputeDiffusionTerm,
            " ", py::arg("rX"), py::arg("pElement") = __null)
        .def("GetUptakeRateForElement",
            (double(AveragedSourceParabolicPde_3::*)(unsigned int)) &AveragedSourceParabolicPde_3::GetUptakeRateForElement,
            " ", py::arg("elementIndex"))
    ;
}
