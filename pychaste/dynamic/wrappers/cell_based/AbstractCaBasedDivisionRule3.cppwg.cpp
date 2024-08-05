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
#include "AbstractCaBasedDivisionRule.hpp"

#include "AbstractCaBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaBasedDivisionRule<3 > AbstractCaBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractCaBasedDivisionRule3_Overrides : public AbstractCaBasedDivisionRule3{
    public:
    using AbstractCaBasedDivisionRule3::AbstractCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCaBasedDivisionRule3,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractCaBasedDivisionRule3,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }
    void OutputCellCaBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCaBasedDivisionRule3,
            OutputCellCaBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaBasedDivisionRule3_class(py::module &m){
py::class_<AbstractCaBasedDivisionRule3 , AbstractCaBasedDivisionRule3_Overrides , boost::shared_ptr<AbstractCaBasedDivisionRule3 >   >(m, "AbstractCaBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "IsRoomToDivide",
            (bool(AbstractCaBasedDivisionRule3::*)(::CellPtr, ::CaBasedCellPopulation<3> &)) &AbstractCaBasedDivisionRule3::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(AbstractCaBasedDivisionRule3::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<3> &)) &AbstractCaBasedDivisionRule3::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCaBasedDivisionRuleInfo",
            (void(AbstractCaBasedDivisionRule3::*)(::out_stream &)) &AbstractCaBasedDivisionRule3::OutputCellCaBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
