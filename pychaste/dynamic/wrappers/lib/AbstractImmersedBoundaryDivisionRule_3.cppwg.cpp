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
#include "AbstractImmersedBoundaryDivisionRule.hpp"

#include "AbstractImmersedBoundaryDivisionRule_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractImmersedBoundaryDivisionRule<3> AbstractImmersedBoundaryDivisionRule_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class AbstractImmersedBoundaryDivisionRule_3_Overrides : public AbstractImmersedBoundaryDivisionRule_3
{
public:
    using AbstractImmersedBoundaryDivisionRule_3::AbstractImmersedBoundaryDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::ImmersedBoundaryCellPopulation<3> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractImmersedBoundaryDivisionRule_3,
            CalculateCellDivisionVector,
            pParentCell,
            rCellPopulation);
    }
    void OutputCellImmersedBoundaryDivisionRuleParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractImmersedBoundaryDivisionRule_3,
            OutputCellImmersedBoundaryDivisionRuleParameters,
            rParamsFile);
    }
};

void register_AbstractImmersedBoundaryDivisionRule_3_class(py::module &m)
{
    py::class_<AbstractImmersedBoundaryDivisionRule_3, AbstractImmersedBoundaryDivisionRule_3_Overrides, boost::shared_ptr<AbstractImmersedBoundaryDivisionRule_3>, Identifiable>(m, "AbstractImmersedBoundaryDivisionRule_3")
        .def(py::init<>())
        .def("CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractImmersedBoundaryDivisionRule_3::*)(::CellPtr, ::ImmersedBoundaryCellPopulation<3> &)) &AbstractImmersedBoundaryDivisionRule_3::CalculateCellDivisionVector,
            " ", py::arg("pParentCell"), py::arg("rCellPopulation"))
        .def("OutputCellImmersedBoundaryDivisionRuleInfo",
            (void(AbstractImmersedBoundaryDivisionRule_3::*)(::out_stream &)) &AbstractImmersedBoundaryDivisionRule_3::OutputCellImmersedBoundaryDivisionRuleInfo,
            " ", py::arg("rParamsFile"))
    ;
}
