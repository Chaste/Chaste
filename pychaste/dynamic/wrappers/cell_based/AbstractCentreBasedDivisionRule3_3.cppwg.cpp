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
#include "AbstractCentreBasedDivisionRule.hpp"

#include "AbstractCentreBasedDivisionRule3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCentreBasedDivisionRule<3,3 > AbstractCentreBasedDivisionRule3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_;

class AbstractCentreBasedDivisionRule3_3_Overrides : public AbstractCentreBasedDivisionRule3_3{
    public:
    using AbstractCentreBasedDivisionRule3_3::AbstractCentreBasedDivisionRule;
    ::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>> CalculateCellDivisionVector(::CellPtr pParentCell, ::AbstractCentreBasedCellPopulation<3, 3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _std_pair_lt_boost_numeric_ublas_c_vector_lt_double_3_gt__boost_numeric_ublas_c_vector_lt_double_3_gt__gt_,
            AbstractCentreBasedDivisionRule3_3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellCentreBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedDivisionRule3_3,
            OutputCellCentreBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCentreBasedDivisionRule3_3_class(py::module &m){
py::class_<AbstractCentreBasedDivisionRule3_3 , AbstractCentreBasedDivisionRule3_3_Overrides , boost::shared_ptr<AbstractCentreBasedDivisionRule3_3 >   >(m, "AbstractCentreBasedDivisionRule3_3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::std::pair<boost::numeric::ublas::c_vector<double, 3>, boost::numeric::ublas::c_vector<double, 3>>(AbstractCentreBasedDivisionRule3_3::*)(::CellPtr, ::AbstractCentreBasedCellPopulation<3, 3> &)) &AbstractCentreBasedDivisionRule3_3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCentreBasedDivisionRuleInfo",
            (void(AbstractCentreBasedDivisionRule3_3::*)(::out_stream &)) &AbstractCentreBasedDivisionRule3_3::OutputCellCentreBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
