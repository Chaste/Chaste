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
#include "LabelDependentBernoulliTrialCellCycleModel.hpp"

#include "LabelDependentBernoulliTrialCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef LabelDependentBernoulliTrialCellCycleModel LabelDependentBernoulliTrialCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class LabelDependentBernoulliTrialCellCycleModel_Overrides : public LabelDependentBernoulliTrialCellCycleModel
{
public:
    using LabelDependentBernoulliTrialCellCycleModel::LabelDependentBernoulliTrialCellCycleModel;
    bool ReadyToDivide() override
    {
        PYBIND11_OVERRIDE(
            bool,
            LabelDependentBernoulliTrialCellCycleModel,
            ReadyToDivide,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override
    {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            LabelDependentBernoulliTrialCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override
    {
        PYBIND11_OVERRIDE(
            double,
            LabelDependentBernoulliTrialCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override
    {
        PYBIND11_OVERRIDE(
            double,
            LabelDependentBernoulliTrialCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            LabelDependentBernoulliTrialCellCycleModel,
            OutputCellCycleModelParameters,
            rParamsFile);
    }
};

void register_LabelDependentBernoulliTrialCellCycleModel_class(py::module &m)
{
    py::class_<LabelDependentBernoulliTrialCellCycleModel, LabelDependentBernoulliTrialCellCycleModel_Overrides, boost::shared_ptr<LabelDependentBernoulliTrialCellCycleModel>, AbstractCellCycleModel>(m, "LabelDependentBernoulliTrialCellCycleModel")
        .def(py::init<>())
        .def("ReadyToDivide",
            (bool(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::ReadyToDivide,
            " ")
        .def("CreateCellCycleModel",
            (::AbstractCellCycleModel *(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::CreateCellCycleModel,
            " ", py::return_value_policy::reference)
        .def("SetDivisionProbability",
            (void(LabelDependentBernoulliTrialCellCycleModel::*)(double)) &LabelDependentBernoulliTrialCellCycleModel::SetDivisionProbability,
            " ", py::arg("divisionProbability"))
        .def("GetLabelledDivisionProbability",
            (double(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::GetLabelledDivisionProbability,
            " ")
        .def("SetLabelledDivisionProbability",
            (void(LabelDependentBernoulliTrialCellCycleModel::*)(double)) &LabelDependentBernoulliTrialCellCycleModel::SetLabelledDivisionProbability,
            " ", py::arg("labelledDivisionProbability"))
        .def("GetDivisionProbability",
            (double(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::GetDivisionProbability,
            " ")
        .def("SetMinimumDivisionAge",
            (void(LabelDependentBernoulliTrialCellCycleModel::*)(double)) &LabelDependentBernoulliTrialCellCycleModel::SetMinimumDivisionAge,
            " ", py::arg("minimumDivisionAge"))
        .def("GetMinimumDivisionAge",
            (double(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::GetMinimumDivisionAge,
            " ")
        .def("GetAverageTransitCellCycleTime",
            (double(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::GetAverageTransitCellCycleTime,
            " ")
        .def("GetAverageStemCellCycleTime",
            (double(LabelDependentBernoulliTrialCellCycleModel::*)()) &LabelDependentBernoulliTrialCellCycleModel::GetAverageStemCellCycleTime,
            " ")
        .def("OutputCellCycleModelParameters",
            (void(LabelDependentBernoulliTrialCellCycleModel::*)(::out_stream &)) &LabelDependentBernoulliTrialCellCycleModel::OutputCellCycleModelParameters,
            " ", py::arg("rParamsFile"))
    ;
}
