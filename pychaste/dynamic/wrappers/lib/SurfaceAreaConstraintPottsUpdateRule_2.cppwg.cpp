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
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "SurfaceAreaConstraintPottsUpdateRule_2.cppwg.hpp"

namespace py = pybind11;
typedef SurfaceAreaConstraintPottsUpdateRule<2> SurfaceAreaConstraintPottsUpdateRule_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SurfaceAreaConstraintPottsUpdateRule_2_Overrides : public SurfaceAreaConstraintPottsUpdateRule_2
{
public:
    using SurfaceAreaConstraintPottsUpdateRule_2::SurfaceAreaConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            double,
            SurfaceAreaConstraintPottsUpdateRule_2,
            EvaluateHamiltonianContribution,
            currentNodeIndex,
            targetNodeIndex,
            rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            SurfaceAreaConstraintPottsUpdateRule_2,
            OutputUpdateRuleParameters,
            rParamsFile);
    }
};

void register_SurfaceAreaConstraintPottsUpdateRule_2_class(py::module &m)
{
    py::class_<SurfaceAreaConstraintPottsUpdateRule_2, SurfaceAreaConstraintPottsUpdateRule_2_Overrides, boost::shared_ptr<SurfaceAreaConstraintPottsUpdateRule_2>, AbstractPottsUpdateRule<2>>(m, "SurfaceAreaConstraintPottsUpdateRule_2")
        .def(py::init<>())
        .def("EvaluateHamiltonianContribution",
            (double(SurfaceAreaConstraintPottsUpdateRule_2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &SurfaceAreaConstraintPottsUpdateRule_2::EvaluateHamiltonianContribution,
            " ", py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation"))
        .def("GetDeformationEnergyParameter",
            (double(SurfaceAreaConstraintPottsUpdateRule_2::*)()) &SurfaceAreaConstraintPottsUpdateRule_2::GetDeformationEnergyParameter,
            " ")
        .def("SetDeformationEnergyParameter",
            (void(SurfaceAreaConstraintPottsUpdateRule_2::*)(double)) &SurfaceAreaConstraintPottsUpdateRule_2::SetDeformationEnergyParameter,
            " ", py::arg("deformationEnergyParameter"))
        .def("GetMatureCellTargetSurfaceArea",
            (double(SurfaceAreaConstraintPottsUpdateRule_2::*)() const) &SurfaceAreaConstraintPottsUpdateRule_2::GetMatureCellTargetSurfaceArea,
            " ")
        .def("SetMatureCellTargetSurfaceArea",
            (void(SurfaceAreaConstraintPottsUpdateRule_2::*)(double)) &SurfaceAreaConstraintPottsUpdateRule_2::SetMatureCellTargetSurfaceArea,
            " ", py::arg("matureCellTargetSurfaceArea"))
        .def("OutputUpdateRuleParameters",
            (void(SurfaceAreaConstraintPottsUpdateRule_2::*)(::out_stream &)) &SurfaceAreaConstraintPottsUpdateRule_2::OutputUpdateRuleParameters,
            " ", py::arg("rParamsFile"))
    ;
}
