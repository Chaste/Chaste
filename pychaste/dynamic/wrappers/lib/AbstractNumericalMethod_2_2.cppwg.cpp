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
#include "AbstractNumericalMethod.hpp"

#include "AbstractNumericalMethod_2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNumericalMethod<2, 2> AbstractNumericalMethod_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractNumericalMethod_2_2_Overrides : public AbstractNumericalMethod_2_2
{
public:
    using AbstractNumericalMethod_2_2::AbstractNumericalMethod;
    void UpdateAllNodePositions(double dt) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractNumericalMethod_2_2,
            UpdateAllNodePositions,
            dt);
    }
    void OutputNumericalMethodParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractNumericalMethod_2_2,
            OutputNumericalMethodParameters,
            rParamsFile);
    }
};

void register_AbstractNumericalMethod_2_2_class(py::module &m)
{
    py::class_<AbstractNumericalMethod_2_2, AbstractNumericalMethod_2_2_Overrides, boost::shared_ptr<AbstractNumericalMethod_2_2>, Identifiable>(m, "AbstractNumericalMethod_2_2")
        .def(py::init<>())
        .def("SetCellPopulation",
            (void(AbstractNumericalMethod_2_2::*)(::AbstractOffLatticeCellPopulation<2> *)) &AbstractNumericalMethod_2_2::SetCellPopulation,
            " ", py::arg("pPopulation"))
        .def("SetForceCollection",
            (void(AbstractNumericalMethod_2_2::*)(::std::vector<boost::shared_ptr<AbstractForce<2, 2>>> *)) &AbstractNumericalMethod_2_2::SetForceCollection,
            " ", py::arg("pForces"))
        .def("SetBoundaryConditions",
            (void(AbstractNumericalMethod_2_2::*)(::std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<2, 2>>> *)) &AbstractNumericalMethod_2_2::SetBoundaryConditions,
            " ", py::arg("pBoundaryConditions"))
        .def("SetUseAdaptiveTimestep",
            (void(AbstractNumericalMethod_2_2::*)(bool)) &AbstractNumericalMethod_2_2::SetUseAdaptiveTimestep,
            " ", py::arg("useAdaptiveTimestep"))
        .def("SetUseUpdateNodeLocation",
            (void(AbstractNumericalMethod_2_2::*)(bool)) &AbstractNumericalMethod_2_2::SetUseUpdateNodeLocation,
            " ", py::arg("useUpdateNodeLocation"))
        .def("GetUseUpdateNodeLocation",
            (bool(AbstractNumericalMethod_2_2::*)()) &AbstractNumericalMethod_2_2::GetUseUpdateNodeLocation,
            " ")
        .def("HasAdaptiveTimestep",
            (bool(AbstractNumericalMethod_2_2::*)()) &AbstractNumericalMethod_2_2::HasAdaptiveTimestep,
            " ")
        .def("UpdateAllNodePositions",
            (void(AbstractNumericalMethod_2_2::*)(double)) &AbstractNumericalMethod_2_2::UpdateAllNodePositions,
            " ", py::arg("dt"))
        .def("OutputNumericalMethodInfo",
            (void(AbstractNumericalMethod_2_2::*)(::out_stream &)) &AbstractNumericalMethod_2_2::OutputNumericalMethodInfo,
            " ", py::arg("rParamsFile"))
        .def("OutputNumericalMethodParameters",
            (void(AbstractNumericalMethod_2_2::*)(::out_stream &)) &AbstractNumericalMethod_2_2::OutputNumericalMethodParameters,
            " ", py::arg("rParamsFile"))
    ;
}
