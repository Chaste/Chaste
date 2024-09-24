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
#include "AbstractCaUpdateRule.hpp"

#include "AbstractCaUpdateRule_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaUpdateRule<3> AbstractCaUpdateRule_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaUpdateRule_3_Overrides : public AbstractCaUpdateRule_3
{
public:
    using AbstractCaUpdateRule_3::AbstractCaUpdateRule;
    double EvaluateProbability(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::CaBasedCellPopulation<3> & rCellPopulation, double dt, double deltaX, ::CellPtr cell) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaUpdateRule_3,
            EvaluateProbability,
            currentNodeIndex,
            targetNodeIndex,
            rCellPopulation,
            dt,
            deltaX,
            cell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCaUpdateRule_3,
            OutputUpdateRuleParameters,
            rParamsFile);
    }
};

void register_AbstractCaUpdateRule_3_class(py::module &m)
{
    py::class_<AbstractCaUpdateRule_3, AbstractCaUpdateRule_3_Overrides, boost::shared_ptr<AbstractCaUpdateRule_3>, AbstractUpdateRule<3>>(m, "AbstractCaUpdateRule_3")
        .def("EvaluateProbability",
            (double(AbstractCaUpdateRule_3::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<3> &, double, double, ::CellPtr)) &AbstractCaUpdateRule_3::EvaluateProbability,
            " ", py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX"), py::arg("cell"))
        .def("OutputUpdateRuleParameters",
            (void(AbstractCaUpdateRule_3::*)(::out_stream &)) &AbstractCaUpdateRule_3::OutputUpdateRuleParameters,
            " ", py::arg("rParamsFile"))
    ;
}
