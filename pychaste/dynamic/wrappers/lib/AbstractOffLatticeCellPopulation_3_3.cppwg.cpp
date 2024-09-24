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
#include "AbstractOffLatticeCellPopulation.hpp"

#include "AbstractOffLatticeCellPopulation_3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOffLatticeCellPopulation<3, 3> AbstractOffLatticeCellPopulation_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractOffLatticeCellPopulation_3_3_Overrides : public AbstractOffLatticeCellPopulation_3_3
{
public:
    using AbstractOffLatticeCellPopulation_3_3::AbstractOffLatticeCellPopulation;
    unsigned int AddNode(::Node<3> * pNewNode) override
    {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractOffLatticeCellPopulation_3_3,
            AddNode,
            pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> & rNewLocation) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation_3_3,
            SetNode,
            nodeIndex,
            rNewLocation);
    }
    void UpdateNodeLocations(double dt) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation_3_3,
            UpdateNodeLocations,
            dt);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 3> & rDisplacement, double dt) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation_3_3,
            CheckForStepSizeException,
            nodeIndex,
            rDisplacement,
            dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractOffLatticeCellPopulation_3_3,
            GetDampingConstant,
            nodeIndex);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation_3_3,
            OutputCellPopulationParameters,
            rParamsFile);
    }
};

void register_AbstractOffLatticeCellPopulation_3_3_class(py::module &m)
{
    py::class_<AbstractOffLatticeCellPopulation_3_3, AbstractOffLatticeCellPopulation_3_3_Overrides, boost::shared_ptr<AbstractOffLatticeCellPopulation_3_3>, AbstractCellPopulation<3>>(m, "AbstractOffLatticeCellPopulation_3_3")
        .def("AddNode",
            (unsigned int(AbstractOffLatticeCellPopulation_3_3::*)(::Node<3> *)) &AbstractOffLatticeCellPopulation_3_3::AddNode,
            " ", py::arg("pNewNode"))
        .def("SetNode",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(unsigned int, ::ChastePoint<3> &)) &AbstractOffLatticeCellPopulation_3_3::SetNode,
            " ", py::arg("nodeIndex"), py::arg("rNewLocation"))
        .def("UpdateNodeLocations",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(double)) &AbstractOffLatticeCellPopulation_3_3::UpdateNodeLocations,
            " ", py::arg("dt"))
        .def("CheckForStepSizeException",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 3> &, double)) &AbstractOffLatticeCellPopulation_3_3::CheckForStepSizeException,
            " ", py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt"))
        .def("GetDampingConstant",
            (double(AbstractOffLatticeCellPopulation_3_3::*)(unsigned int)) &AbstractOffLatticeCellPopulation_3_3::GetDampingConstant,
            " ", py::arg("nodeIndex"))
        .def("SetDampingConstantNormal",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(double)) &AbstractOffLatticeCellPopulation_3_3::SetDampingConstantNormal,
            " ", py::arg("dampingConstantNormal"))
        .def("SetDampingConstantMutant",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(double)) &AbstractOffLatticeCellPopulation_3_3::SetDampingConstantMutant,
            " ", py::arg("dampingConstantMutant"))
        .def("SetAbsoluteMovementThreshold",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(double)) &AbstractOffLatticeCellPopulation_3_3::SetAbsoluteMovementThreshold,
            " ", py::arg("absoluteMovementThreshold"))
        .def("GetAbsoluteMovementThreshold",
            (double(AbstractOffLatticeCellPopulation_3_3::*)()) &AbstractOffLatticeCellPopulation_3_3::GetAbsoluteMovementThreshold,
            " ")
        .def("GetDampingConstantNormal",
            (double(AbstractOffLatticeCellPopulation_3_3::*)()) &AbstractOffLatticeCellPopulation_3_3::GetDampingConstantNormal,
            " ")
        .def("GetDampingConstantMutant",
            (double(AbstractOffLatticeCellPopulation_3_3::*)()) &AbstractOffLatticeCellPopulation_3_3::GetDampingConstantMutant,
            " ")
        .def("OutputCellPopulationParameters",
            (void(AbstractOffLatticeCellPopulation_3_3::*)(::out_stream &)) &AbstractOffLatticeCellPopulation_3_3::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
    ;
}
