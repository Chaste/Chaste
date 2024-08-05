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
#include "VonMisesVertexBasedDivisionRule.hpp"

#include "VonMisesVertexBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef VonMisesVertexBasedDivisionRule<2 > VonMisesVertexBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class VonMisesVertexBasedDivisionRule2_Overrides : public VonMisesVertexBasedDivisionRule2{
    public:
    using VonMisesVertexBasedDivisionRule2::VonMisesVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            VonMisesVertexBasedDivisionRule2,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }

};
void register_VonMisesVertexBasedDivisionRule2_class(py::module &m){
py::class_<VonMisesVertexBasedDivisionRule2 , VonMisesVertexBasedDivisionRule2_Overrides , boost::shared_ptr<VonMisesVertexBasedDivisionRule2 >  , AbstractVertexBasedDivisionRule<2>  >(m, "VonMisesVertexBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "GetMeanParameter",
            (double(VonMisesVertexBasedDivisionRule2::*)()) &VonMisesVertexBasedDivisionRule2::GetMeanParameter,
            " "  )
        .def(
            "GetConcentrationParameter",
            (double(VonMisesVertexBasedDivisionRule2::*)()) &VonMisesVertexBasedDivisionRule2::GetConcentrationParameter,
            " "  )
        .def(
            "SetMeanParameter",
            (void(VonMisesVertexBasedDivisionRule2::*)(double)) &VonMisesVertexBasedDivisionRule2::SetMeanParameter,
            " " , py::arg("meanParameter") )
        .def(
            "SetConcentrationParameter",
            (void(VonMisesVertexBasedDivisionRule2::*)(double)) &VonMisesVertexBasedDivisionRule2::SetConcentrationParameter,
            " " , py::arg("concentrationParameter") )
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 2>(VonMisesVertexBasedDivisionRule2::*)(::CellPtr, ::VertexBasedCellPopulation<2> &)) &VonMisesVertexBasedDivisionRule2::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
