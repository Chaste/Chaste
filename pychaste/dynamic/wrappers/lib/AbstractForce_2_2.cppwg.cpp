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
#include "AbstractForce.hpp"

#include "AbstractForce_2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractForce<2, 2> AbstractForce_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractForce_2_2_Overrides : public AbstractForce_2_2
{
public:
    using AbstractForce_2_2::AbstractForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractForce_2_2,
            AddForceContribution,
            rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractForce_2_2,
            OutputForceParameters,
            rParamsFile);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractForce_2_2,
            WriteDataToVisualizerSetupFile,
            pVizSetupFile);
    }
};

void register_AbstractForce_2_2_class(py::module &m)
{
    py::class_<AbstractForce_2_2, AbstractForce_2_2_Overrides, boost::shared_ptr<AbstractForce_2_2>, Identifiable>(m, "AbstractForce_2_2")
        .def(py::init<>())
        .def("AddForceContribution",
            (void(AbstractForce_2_2::*)(::AbstractCellPopulation<2> &)) &AbstractForce_2_2::AddForceContribution,
            " ", py::arg("rCellPopulation"))
        .def("OutputForceInfo",
            (void(AbstractForce_2_2::*)(::out_stream &)) &AbstractForce_2_2::OutputForceInfo,
            " ", py::arg("rParamsFile"))
        .def("OutputForceParameters",
            (void(AbstractForce_2_2::*)(::out_stream &)) &AbstractForce_2_2::OutputForceParameters,
            " ", py::arg("rParamsFile"))
        .def("WriteDataToVisualizerSetupFile",
            (void(AbstractForce_2_2::*)(::out_stream &)) &AbstractForce_2_2::WriteDataToVisualizerSetupFile,
            " ", py::arg("pVizSetupFile"))
    ;
}
