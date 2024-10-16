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
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "CellBasedEllipticPdeSolver.hpp"

#include "CellBasedEllipticPdeSolver_3.cppwg.hpp"

namespace py = pybind11;
typedef CellBasedEllipticPdeSolver<3> CellBasedEllipticPdeSolver_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 4> _boost_numeric_ublas_c_vector_lt_double_4_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 4, 4> _boost_numeric_ublas_c_matrix_lt_double_4_4_gt_;

class CellBasedEllipticPdeSolver_3_Overrides : public CellBasedEllipticPdeSolver_3
{
public:
    using CellBasedEllipticPdeSolver_3::CellBasedEllipticPdeSolver;
    ::boost::numeric::ublas::c_vector<double, 4> ComputeVectorTerm(::boost::numeric::ublas::c_vector<double, 4> & rPhi, ::boost::numeric::ublas::c_matrix<double, 3, 4> & rGradPhi, ::ChastePoint<3> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 3> & rGradU, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_4_gt_,
            CellBasedEllipticPdeSolver_3,
            ComputeVectorTerm,
            rPhi,
            rGradPhi,
            rX,
            rU,
            rGradU,
            pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 4, 4> ComputeMatrixTerm(::boost::numeric::ublas::c_vector<double, 4> & rPhi, ::boost::numeric::ublas::c_matrix<double, 3, 4> & rGradPhi, ::ChastePoint<3> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 3> & rGradU, ::Element<3, 3> * pElement) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_4_4_gt_,
            CellBasedEllipticPdeSolver_3,
            ComputeMatrixTerm,
            rPhi,
            rGradPhi,
            rX,
            rU,
            rGradU,
            pElement);
    }
    void ResetInterpolatedQuantities() override
    {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver_3,
            ResetInterpolatedQuantities,
            );
    }
    void IncrementInterpolatedQuantities(double phiI, ::Node<3> const * pNode) override
    {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver_3,
            IncrementInterpolatedQuantities,
            phiI,
            pNode);
    }
    void InitialiseForSolve(::Vec initialSolution) override
    {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver_3,
            InitialiseForSolve,
            initialSolution);
    }
};

void register_CellBasedEllipticPdeSolver_3_class(py::module &m)
{
    py::class_<CellBasedEllipticPdeSolver_3, CellBasedEllipticPdeSolver_3_Overrides, boost::shared_ptr<CellBasedEllipticPdeSolver_3>>(m, "CellBasedEllipticPdeSolver_3")
        .def(py::init<::TetrahedralMesh<3, 3> *, ::AbstractLinearEllipticPde<3, 3> *, ::BoundaryConditionsContainer<3, 3, 1> *>(), py::arg("pMesh"), py::arg("pPde"), py::arg("pBoundaryConditions"))
    ;
}
