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
#include "PythonUblasObjectConverters.hpp"
#include "ApoptoticCellKiller.hpp"

#include "ApoptoticCellKiller_2.cppwg.hpp"

namespace py = pybind11;
typedef ApoptoticCellKiller<2> ApoptoticCellKiller_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ApoptoticCellKiller_2_Overrides : public ApoptoticCellKiller_2
{
public:
    using ApoptoticCellKiller_2::ApoptoticCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override
    {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller_2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller_2,
            OutputCellKillerParameters,
            rParamsFile);
    }
};

void register_ApoptoticCellKiller_2_class(py::module &m)
{
    py::class_<ApoptoticCellKiller_2, ApoptoticCellKiller_2_Overrides, boost::shared_ptr<ApoptoticCellKiller_2>, AbstractCellKiller<2>>(m, "ApoptoticCellKiller_2")
        .def(py::init<::AbstractCellPopulation<2> *>(), py::arg("pCellPopulation"))
        .def("CheckAndLabelSingleCellForApoptosis",
            (void(ApoptoticCellKiller_2::*)(::CellPtr)) &ApoptoticCellKiller_2::CheckAndLabelSingleCellForApoptosis,
            " ", py::arg("pCell"))
        .def("CheckAndLabelCellsForApoptosisOrDeath",
            (void(ApoptoticCellKiller_2::*)()) &ApoptoticCellKiller_2::CheckAndLabelCellsForApoptosisOrDeath,
            " ")
        .def("OutputCellKillerParameters",
            (void(ApoptoticCellKiller_2::*)(::out_stream &)) &ApoptoticCellKiller_2::OutputCellKillerParameters,
            " ", py::arg("rParamsFile"))
    ;
}
