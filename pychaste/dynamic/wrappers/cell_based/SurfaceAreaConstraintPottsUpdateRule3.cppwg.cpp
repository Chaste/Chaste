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
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "SurfaceAreaConstraintPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef SurfaceAreaConstraintPottsUpdateRule<3 > SurfaceAreaConstraintPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SurfaceAreaConstraintPottsUpdateRule3_Overrides : public SurfaceAreaConstraintPottsUpdateRule3{
    public:
    using SurfaceAreaConstraintPottsUpdateRule3::SurfaceAreaConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            SurfaceAreaConstraintPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SurfaceAreaConstraintPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_SurfaceAreaConstraintPottsUpdateRule3_class(py::module &m){
py::class_<SurfaceAreaConstraintPottsUpdateRule3 , SurfaceAreaConstraintPottsUpdateRule3_Overrides , boost::shared_ptr<SurfaceAreaConstraintPottsUpdateRule3 >  , AbstractPottsUpdateRule<3>  >(m, "SurfaceAreaConstraintPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &SurfaceAreaConstraintPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetDeformationEnergyParameter",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)()) &SurfaceAreaConstraintPottsUpdateRule3::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(double)) &SurfaceAreaConstraintPottsUpdateRule3::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "GetMatureCellTargetSurfaceArea",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)() const ) &SurfaceAreaConstraintPottsUpdateRule3::GetMatureCellTargetSurfaceArea,
            " "  )
        .def(
            "SetMatureCellTargetSurfaceArea",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(double)) &SurfaceAreaConstraintPottsUpdateRule3::SetMatureCellTargetSurfaceArea,
            " " , py::arg("matureCellTargetSurfaceArea") )
        .def(
            "OutputUpdateRuleParameters",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(::out_stream &)) &SurfaceAreaConstraintPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
