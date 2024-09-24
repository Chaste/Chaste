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
#include "AbstractCellBasedSimulation.hpp"

#include "AbstractCellBasedSimulation_3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedSimulation<3, 3> AbstractCellBasedSimulation_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractCellBasedSimulation_3_3_Overrides : public AbstractCellBasedSimulation_3_3
{
public:
    using AbstractCellBasedSimulation_3_3::AbstractCellBasedSimulation;
    void OutputSimulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation_3_3,
            OutputSimulationParameters,
            rParamsFile);
    }
    void WriteVisualizerSetupFile() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation_3_3,
            WriteVisualizerSetupFile,
            );
    }
    unsigned int DoCellBirth() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractCellBasedSimulation_3_3,
            DoCellBirth,
            );
    }
    void SetupSolve() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation_3_3,
            SetupSolve,
            );
    }
    bool StoppingEventHasOccurred() override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellBasedSimulation_3_3,
            StoppingEventHasOccurred,
            );
    }
    void UpdateCellPopulation() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation_3_3,
            UpdateCellPopulation,
            );
    }
    void UpdateCellLocationsAndTopology() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation_3_3,
            UpdateCellLocationsAndTopology,
            );
    }
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation_3_3,
            OutputAdditionalSimulationSetup,
            rParamsFile);
    }
};

void register_AbstractCellBasedSimulation_3_3_class(py::module &m)
{
    py::class_<AbstractCellBasedSimulation_3_3, AbstractCellBasedSimulation_3_3_Overrides, boost::shared_ptr<AbstractCellBasedSimulation_3_3>, Identifiable>(m, "AbstractCellBasedSimulation_3_3")
        .def(py::init<::AbstractCellPopulation<3> &, bool, bool>(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def("GetNodeLocation",
            (::std::vector<double>(AbstractCellBasedSimulation_3_3::*)(unsigned int const &)) &AbstractCellBasedSimulation_3_3::GetNodeLocation,
            " ", py::arg("rNodeIndex"))
        .def("GetDt",
            (double(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetDt,
            " ")
        .def("GetNumBirths",
            (unsigned int(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetNumBirths,
            " ")
        .def("GetNumDeaths",
            (unsigned int(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetNumDeaths,
            " ")
        .def("GetOutputDirectory",
            (::std::string(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetOutputDirectory,
            " ")
        .def("SetDt",
            (void(AbstractCellBasedSimulation_3_3::*)(double)) &AbstractCellBasedSimulation_3_3::SetDt,
            " ", py::arg("dt"))
        .def("SetEndTime",
            (void(AbstractCellBasedSimulation_3_3::*)(double)) &AbstractCellBasedSimulation_3_3::SetEndTime,
            " ", py::arg("endTime"))
        .def("SetOutputDirectory",
            (void(AbstractCellBasedSimulation_3_3::*)(::std::string)) &AbstractCellBasedSimulation_3_3::SetOutputDirectory,
            " ", py::arg("outputDirectory"))
        .def("SetSamplingTimestepMultiple",
            (void(AbstractCellBasedSimulation_3_3::*)(unsigned int)) &AbstractCellBasedSimulation_3_3::SetSamplingTimestepMultiple,
            " ", py::arg("samplingTimestepMultiple"))
        .def("SetUpdatingTimestepMultiple",
            (void(AbstractCellBasedSimulation_3_3::*)(unsigned int)) &AbstractCellBasedSimulation_3_3::SetUpdatingTimestepMultiple,
            " ", py::arg("updatingTimestepMultiple"))
        .def("SetNoBirth",
            (void(AbstractCellBasedSimulation_3_3::*)(bool)) &AbstractCellBasedSimulation_3_3::SetNoBirth,
            " ", py::arg("noBirth"))
        .def("SetUpdateCellPopulationRule",
            (void(AbstractCellBasedSimulation_3_3::*)(bool)) &AbstractCellBasedSimulation_3_3::SetUpdateCellPopulationRule,
            " ", py::arg("updateCellPopulation"))
        .def("GetUpdateCellPopulationRule",
            (bool(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetUpdateCellPopulationRule,
            " ")
        .def("AddCellKiller",
            (void(AbstractCellBasedSimulation_3_3::*)(::boost::shared_ptr<AbstractCellKiller<3>>)) &AbstractCellBasedSimulation_3_3::AddCellKiller,
            " ", py::arg("pCellKiller"))
        .def("RemoveAllCellKillers",
            (void(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::RemoveAllCellKillers,
            " ")
        .def("AddSimulationModifier",
            (void(AbstractCellBasedSimulation_3_3::*)(::boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>)) &AbstractCellBasedSimulation_3_3::AddSimulationModifier,
            " ", py::arg("pSimulationModifier"))
        .def("AddTopologyUpdateSimulationModifier",
            (void(AbstractCellBasedSimulation_3_3::*)(::boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>)) &AbstractCellBasedSimulation_3_3::AddTopologyUpdateSimulationModifier,
            " ", py::arg("pSimulationModifier"))
        .def("GetTopologyUpdateSimulationModifiers",
            (::std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>> *(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetTopologyUpdateSimulationModifiers,
            " ", py::return_value_policy::reference)
        .def("Solve",
            (void(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::Solve,
            " ")
        .def("rGetCellPopulation",
            (::AbstractCellPopulation<3> &(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::rGetCellPopulation,
            " ", py::return_value_policy::reference_internal)
        .def("rGetCellPopulation",
            (::AbstractCellPopulation<3> const &(AbstractCellBasedSimulation_3_3::*)() const) &AbstractCellBasedSimulation_3_3::rGetCellPopulation,
            " ", py::return_value_policy::reference_internal)
        .def("GetOutputDivisionLocations",
            (bool(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetOutputDivisionLocations,
            " ")
        .def("SetOutputDivisionLocations",
            (void(AbstractCellBasedSimulation_3_3::*)(bool)) &AbstractCellBasedSimulation_3_3::SetOutputDivisionLocations,
            " ", py::arg("outputDivisionLocations"))
        .def("GetOutputCellVelocities",
            (bool(AbstractCellBasedSimulation_3_3::*)()) &AbstractCellBasedSimulation_3_3::GetOutputCellVelocities,
            " ")
        .def("SetOutputCellVelocities",
            (void(AbstractCellBasedSimulation_3_3::*)(bool)) &AbstractCellBasedSimulation_3_3::SetOutputCellVelocities,
            " ", py::arg("outputCellVelocities"))
        .def("OutputSimulationParameters",
            (void(AbstractCellBasedSimulation_3_3::*)(::out_stream &)) &AbstractCellBasedSimulation_3_3::OutputSimulationParameters,
            " ", py::arg("rParamsFile"))
    ;
}
