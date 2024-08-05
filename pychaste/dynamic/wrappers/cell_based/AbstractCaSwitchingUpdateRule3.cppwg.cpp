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
#include "AbstractCaSwitchingUpdateRule.hpp"

#include "AbstractCaSwitchingUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaSwitchingUpdateRule<3 > AbstractCaSwitchingUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaSwitchingUpdateRule3_Overrides : public AbstractCaSwitchingUpdateRule3{
    public:
    using AbstractCaSwitchingUpdateRule3::AbstractCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<3> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaSwitchingUpdateRule3,
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
            AbstractCaSwitchingUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaSwitchingUpdateRule3_class(py::module &m){
py::class_<AbstractCaSwitchingUpdateRule3 , AbstractCaSwitchingUpdateRule3_Overrides , boost::shared_ptr<AbstractCaSwitchingUpdateRule3 >  , AbstractUpdateRule<3>  >(m, "AbstractCaSwitchingUpdateRule3")
        .def(
            "EvaluateSwitchingProbability",
            (double(AbstractCaSwitchingUpdateRule3::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<3> &, double, double)) &AbstractCaSwitchingUpdateRule3::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractCaSwitchingUpdateRule3::*)(::out_stream &)) &AbstractCaSwitchingUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
