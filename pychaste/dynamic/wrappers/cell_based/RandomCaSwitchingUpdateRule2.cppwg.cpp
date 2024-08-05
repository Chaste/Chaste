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
#include "RandomCaSwitchingUpdateRule.hpp"

#include "RandomCaSwitchingUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef RandomCaSwitchingUpdateRule<2 > RandomCaSwitchingUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RandomCaSwitchingUpdateRule2_Overrides : public RandomCaSwitchingUpdateRule2{
    public:
    using RandomCaSwitchingUpdateRule2::RandomCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<2> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE(
            double,
            RandomCaSwitchingUpdateRule2,
            EvaluateSwitchingProbability,
                    currentNodeIndex,
        neighbourNodeIndex,
        rCellPopulation,
        dt,
        deltaX);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            RandomCaSwitchingUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_RandomCaSwitchingUpdateRule2_class(py::module &m){
py::class_<RandomCaSwitchingUpdateRule2 , RandomCaSwitchingUpdateRule2_Overrides , boost::shared_ptr<RandomCaSwitchingUpdateRule2 >  , AbstractCaSwitchingUpdateRule<2>  >(m, "RandomCaSwitchingUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateSwitchingProbability",
            (double(RandomCaSwitchingUpdateRule2::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<2> &, double, double)) &RandomCaSwitchingUpdateRule2::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "GetSwitchingParameter",
            (double(RandomCaSwitchingUpdateRule2::*)()) &RandomCaSwitchingUpdateRule2::GetSwitchingParameter,
            " "  )
        .def(
            "SetSwitchingParameter",
            (void(RandomCaSwitchingUpdateRule2::*)(double)) &RandomCaSwitchingUpdateRule2::SetSwitchingParameter,
            " " , py::arg("switchingParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(RandomCaSwitchingUpdateRule2::*)(::out_stream &)) &RandomCaSwitchingUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
