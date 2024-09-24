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

#include "AbstractOffLatticeCellPopulation_2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOffLatticeCellPopulation<2, 2> AbstractOffLatticeCellPopulation_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractOffLatticeCellPopulation_2_2_Overrides : public AbstractOffLatticeCellPopulation_2_2
{
public:
    using AbstractOffLatticeCellPopulation_2_2::AbstractOffLatticeCellPopulation;
    unsigned int AddNode(::Node<2> * pNewNode) override
    {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractOffLatticeCellPopulation_2_2,
            AddNode,
            pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> & rNewLocation) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation_2_2,
            SetNode,
            nodeIndex,
            rNewLocation);
    }
    void UpdateNodeLocations(double dt) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation_2_2,
            UpdateNodeLocations,
            dt);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 2> & rDisplacement, double dt) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation_2_2,
            CheckForStepSizeException,
            nodeIndex,
            rDisplacement,
            dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override
    {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractOffLatticeCellPopulation_2_2,
            GetDampingConstant,
            nodeIndex);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation_2_2,
            OutputCellPopulationParameters,
            rParamsFile);
    }
};

void register_AbstractOffLatticeCellPopulation_2_2_class(py::module &m)
{
    py::class_<AbstractOffLatticeCellPopulation_2_2, AbstractOffLatticeCellPopulation_2_2_Overrides, boost::shared_ptr<AbstractOffLatticeCellPopulation_2_2>, AbstractCellPopulation<2>>(m, "AbstractOffLatticeCellPopulation_2_2")
        .def("AddNode",
            (unsigned int(AbstractOffLatticeCellPopulation_2_2::*)(::Node<2> *)) &AbstractOffLatticeCellPopulation_2_2::AddNode,
            " ", py::arg("pNewNode"))
        .def("SetNode",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(unsigned int, ::ChastePoint<2> &)) &AbstractOffLatticeCellPopulation_2_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("rNewLocation"))
        .def("UpdateNodeLocations",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(double)) &AbstractOffLatticeCellPopulation_2_2::UpdateNodeLocations,
            " ", py::arg("dt"))
        .def("CheckForStepSizeException",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2> &, double)) &AbstractOffLatticeCellPopulation_2_2::CheckForStepSizeException,
            " ", py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt"))
        .def("GetDampingConstant",
            (double(AbstractOffLatticeCellPopulation_2_2::*)(unsigned int)) &AbstractOffLatticeCellPopulation_2_2::GetDampingConstant,
            " ", py::arg("nodeIndex"))
        .def("SetDampingConstantNormal",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(double)) &AbstractOffLatticeCellPopulation_2_2::SetDampingConstantNormal,
            " ", py::arg("dampingConstantNormal"))
        .def("SetDampingConstantMutant",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(double)) &AbstractOffLatticeCellPopulation_2_2::SetDampingConstantMutant,
            " ", py::arg("dampingConstantMutant"))
        .def("SetAbsoluteMovementThreshold",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(double)) &AbstractOffLatticeCellPopulation_2_2::SetAbsoluteMovementThreshold,
            " ", py::arg("absoluteMovementThreshold"))
        .def("GetAbsoluteMovementThreshold",
            (double(AbstractOffLatticeCellPopulation_2_2::*)()) &AbstractOffLatticeCellPopulation_2_2::GetAbsoluteMovementThreshold,
            " ")
        .def("GetDampingConstantNormal",
            (double(AbstractOffLatticeCellPopulation_2_2::*)()) &AbstractOffLatticeCellPopulation_2_2::GetDampingConstantNormal,
            " ")
        .def("GetDampingConstantMutant",
            (double(AbstractOffLatticeCellPopulation_2_2::*)()) &AbstractOffLatticeCellPopulation_2_2::GetDampingConstantMutant,
            " ")
        .def("OutputCellPopulationParameters",
            (void(AbstractOffLatticeCellPopulation_2_2::*)(::out_stream &)) &AbstractOffLatticeCellPopulation_2_2::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
    ;
}
