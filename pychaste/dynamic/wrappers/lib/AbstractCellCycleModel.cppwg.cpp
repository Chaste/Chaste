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
#include "AbstractCellCycleModel.hpp"

#include "AbstractCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellCycleModel AbstractCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class AbstractCellCycleModel_Overrides : public AbstractCellCycleModel
{
public:
    using AbstractCellCycleModel::AbstractCellCycleModel;
    void Initialise() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            Initialise,
            );
    }
    void InitialiseDaughterCell() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            InitialiseDaughterCell,
            );
    }
    void SetBirthTime(double birthTime) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            SetBirthTime,
            birthTime);
    }
    bool ReadyToDivide() override
    {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellCycleModel,
            ReadyToDivide,
            );
    }
    void ResetForDivision() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            ResetForDivision,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override
    {
        PYBIND11_OVERRIDE_PURE(
            _AbstractCellCycleModelPtr,
            AbstractCellCycleModel,
            CreateCellCycleModel,
            );
    }
    bool CanCellTerminallyDifferentiate() override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellCycleModel,
            CanCellTerminallyDifferentiate,
            );
    }
    double GetAverageTransitCellCycleTime() override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModel,
            OutputCellCycleModelParameters,
            rParamsFile);
    }
};

void register_AbstractCellCycleModel_class(py::module &m)
{
    py::class_<AbstractCellCycleModel, AbstractCellCycleModel_Overrides, boost::shared_ptr<AbstractCellCycleModel>, Identifiable>(m, "AbstractCellCycleModel")
        .def(py::init<>())
        .def("SetCell",
            (void(AbstractCellCycleModel::*)(::CellPtr)) &AbstractCellCycleModel::SetCell,
            " ", py::arg("pCell"))
        .def("Initialise",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::Initialise,
            " ")
        .def("InitialiseDaughterCell",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::InitialiseDaughterCell,
            " ")
        .def("GetCell",
            (::CellPtr(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetCell,
            " ")
        .def("SetBirthTime",
            (void(AbstractCellCycleModel::*)(double)) &AbstractCellCycleModel::SetBirthTime,
            " ", py::arg("birthTime"))
        .def("SetDimension",
            (void(AbstractCellCycleModel::*)(unsigned int)) &AbstractCellCycleModel::SetDimension,
            " ", py::arg("dimension"))
        .def("GetDimension",
            (unsigned int(AbstractCellCycleModel::*)() const) &AbstractCellCycleModel::GetDimension,
            " ")
        .def("GetBirthTime",
            (double(AbstractCellCycleModel::*)() const) &AbstractCellCycleModel::GetBirthTime,
            " ")
        .def("GetAge",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAge,
            " ")
        .def("ReadyToDivide",
            (bool(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::ReadyToDivide,
            " ")
        .def("ResetForDivision",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::ResetForDivision,
            " ")
        .def("CreateCellCycleModel",
            (::AbstractCellCycleModel *(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::CreateCellCycleModel,
            " ", py::return_value_policy::reference)
        .def("CanCellTerminallyDifferentiate",
            (bool(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::CanCellTerminallyDifferentiate,
            " ")
        .def("GetAverageTransitCellCycleTime",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAverageTransitCellCycleTime,
            " ")
        .def("GetAverageStemCellCycleTime",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAverageStemCellCycleTime,
            " ")
        .def("OutputCellCycleModelInfo",
            (void(AbstractCellCycleModel::*)(::out_stream &)) &AbstractCellCycleModel::OutputCellCycleModelInfo,
            " ", py::arg("rParamsFile"))
        .def("OutputCellCycleModelParameters",
            (void(AbstractCellCycleModel::*)(::out_stream &)) &AbstractCellCycleModel::OutputCellCycleModelParameters,
            " ", py::arg("rParamsFile"))
    ;
}
