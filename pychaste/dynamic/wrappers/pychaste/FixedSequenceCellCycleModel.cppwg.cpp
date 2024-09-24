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
#include "FixedSequenceCellCycleModel.hpp"

#include "FixedSequenceCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef FixedSequenceCellCycleModel FixedSequenceCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class FixedSequenceCellCycleModel_Overrides : public FixedSequenceCellCycleModel
{
public:
    using FixedSequenceCellCycleModel::FixedSequenceCellCycleModel;
    ::AbstractCellCycleModel * CreateCellCycleModel() override
    {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            FixedSequenceCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void SetRate(double rate) override
    {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetRate,
            rate);
    }
    double GetRate() override
    {
        PYBIND11_OVERRIDE(
            double,
            FixedSequenceCellCycleModel,
            GetRate,
            );
    }
    void SetStemCellG1Duration(double stemCellG1Duration) override
    {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetStemCellG1Duration,
            stemCellG1Duration);
    }
    void SetTransitCellG1Duration(double transitCellG1Duration) override
    {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetTransitCellG1Duration,
            transitCellG1Duration);
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            OutputCellCycleModelParameters,
            rParamsFile);
    }
    void SetG1Duration() override
    {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetG1Duration,
            );
    }
};

void register_FixedSequenceCellCycleModel_class(py::module &m)
{
    py::class_<FixedSequenceCellCycleModel, FixedSequenceCellCycleModel_Overrides, boost::shared_ptr<FixedSequenceCellCycleModel>, ExponentialG1GenerationalCellCycleModel>(m, "FixedSequenceCellCycleModel")
        .def(py::init<>())
        .def("CreateCellCycleModel",
            (::AbstractCellCycleModel *(FixedSequenceCellCycleModel::*)()) &FixedSequenceCellCycleModel::CreateCellCycleModel,
            " ", py::return_value_policy::reference)
        .def("SetRate",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetRate,
            " ", py::arg("rate"))
        .def("GetRate",
            (double(FixedSequenceCellCycleModel::*)()) &FixedSequenceCellCycleModel::GetRate,
            " ")
        .def("SetStemCellG1Duration",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetStemCellG1Duration,
            " ", py::arg("stemCellG1Duration"))
        .def("SetTransitCellG1Duration",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetTransitCellG1Duration,
            " ", py::arg("transitCellG1Duration"))
        .def("OutputCellCycleModelParameters",
            (void(FixedSequenceCellCycleModel::*)(::out_stream &)) &FixedSequenceCellCycleModel::OutputCellCycleModelParameters,
            " ", py::arg("rParamsFile"))
    ;
}
