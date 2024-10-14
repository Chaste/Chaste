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
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"

#include "AbstractBoxDomainPdeModifier_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractBoxDomainPdeModifier<2> AbstractBoxDomainPdeModifier_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractBoxDomainPdeModifier_2_Overrides : public AbstractBoxDomainPdeModifier_2
{
public:
    using AbstractBoxDomainPdeModifier_2::AbstractBoxDomainPdeModifier;
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier_2,
            SetupSolve,
            rCellPopulation,
            outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier_2,
            OutputSimulationModifierParameters,
            rParamsFile);
    }
};

void register_AbstractBoxDomainPdeModifier_2_class(py::module &m)
{
    py::class_<AbstractBoxDomainPdeModifier_2, AbstractBoxDomainPdeModifier_2_Overrides, boost::shared_ptr<AbstractBoxDomainPdeModifier_2>, AbstractPdeModifier<2>>(m, "AbstractBoxDomainPdeModifier_2")
        .def("GetStepSize",
            (double(AbstractBoxDomainPdeModifier_2::*)()) &AbstractBoxDomainPdeModifier_2::GetStepSize,
            " ")
        .def("SetBcsOnBoxBoundary",
            (void(AbstractBoxDomainPdeModifier_2::*)(bool)) &AbstractBoxDomainPdeModifier_2::SetBcsOnBoxBoundary,
            " ", py::arg("setBcsOnBoxBoundary"))
        .def("AreBcsSetOnBoxBoundary",
            (bool(AbstractBoxDomainPdeModifier_2::*)()) &AbstractBoxDomainPdeModifier_2::AreBcsSetOnBoxBoundary,
            " ")
        .def("SetupSolve",
            (void(AbstractBoxDomainPdeModifier_2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractBoxDomainPdeModifier_2::SetupSolve,
            " ", py::arg("rCellPopulation"), py::arg("outputDirectory"))
        .def("GenerateFeMesh",
            (void(AbstractBoxDomainPdeModifier_2::*)(::boost::shared_ptr<ChasteCuboid<2>>, double)) &AbstractBoxDomainPdeModifier_2::GenerateFeMesh,
            " ", py::arg("pMeshCuboid"), py::arg("stepSize"))
        .def("UpdateCellData",
            (void(AbstractBoxDomainPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier_2::UpdateCellData,
            " ", py::arg("rCellPopulation"))
        .def("InitialiseCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier_2::InitialiseCellPdeElementMap,
            " ", py::arg("rCellPopulation"))
        .def("UpdateCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier_2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier_2::UpdateCellPdeElementMap,
            " ", py::arg("rCellPopulation"))
        .def("OutputSimulationModifierParameters",
            (void(AbstractBoxDomainPdeModifier_2::*)(::out_stream &)) &AbstractBoxDomainPdeModifier_2::OutputSimulationModifierParameters,
            " ", py::arg("rParamsFile"))
    ;
}
