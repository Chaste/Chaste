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
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"

#include "ParabolicBoxDomainPdeModifier_3.cppwg.hpp"

namespace py = pybind11;
typedef ParabolicBoxDomainPdeModifier<3> ParabolicBoxDomainPdeModifier_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::shared_ptr<BoundaryConditionsContainer<3, 3, 1>> _std_shared_ptr_lt_BoundaryConditionsContainer_lt_3_3_1_gt__gt_;

class ParabolicBoxDomainPdeModifier_3_Overrides : public ParabolicBoxDomainPdeModifier_3
{
public:
    using ParabolicBoxDomainPdeModifier_3::ParabolicBoxDomainPdeModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override
    {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier_3,
            UpdateAtEndOfTimeStep,
            rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier_3,
            SetupSolve,
            rCellPopulation,
            outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier_3,
            OutputSimulationModifierParameters,
            rParamsFile);
    }
};

void register_ParabolicBoxDomainPdeModifier_3_class(py::module &m)
{
    py::class_<ParabolicBoxDomainPdeModifier_3, ParabolicBoxDomainPdeModifier_3_Overrides, boost::shared_ptr<ParabolicBoxDomainPdeModifier_3>, AbstractBoxDomainPdeModifier<3>>(m, "ParabolicBoxDomainPdeModifier_3")
        .def(py::init<::boost::shared_ptr<AbstractLinearPde<3, 3>>, ::boost::shared_ptr<AbstractBoundaryCondition<3>>, bool, ::boost::shared_ptr<ChasteCuboid<3>>, double, ::Vec>(), py::arg("pPde") = boost::shared_ptr<AbstractLinearPde<3, 3>>(), py::arg("pBoundaryCondition") = boost::shared_ptr<AbstractBoundaryCondition<3>>(), py::arg("isNeumannBoundaryCondition") = true, py::arg("pMeshCuboid") = boost::shared_ptr<ChasteCuboid<3>>(), py::arg("stepSize") = 1., py::arg("solution") = nullptr)
        .def("UpdateAtEndOfTimeStep",
            (void(ParabolicBoxDomainPdeModifier_3::*)(::AbstractCellPopulation<3> &)) &ParabolicBoxDomainPdeModifier_3::UpdateAtEndOfTimeStep,
            " ", py::arg("rCellPopulation"))
        .def("SetupSolve",
            (void(ParabolicBoxDomainPdeModifier_3::*)(::AbstractCellPopulation<3> &, ::std::string)) &ParabolicBoxDomainPdeModifier_3::SetupSolve,
            " ", py::arg("rCellPopulation"), py::arg("outputDirectory"))
        .def("SetupInitialSolutionVector",
            (void(ParabolicBoxDomainPdeModifier_3::*)(::AbstractCellPopulation<3> &)) &ParabolicBoxDomainPdeModifier_3::SetupInitialSolutionVector,
            " ", py::arg("rCellPopulation"))
        .def("OutputSimulationModifierParameters",
            (void(ParabolicBoxDomainPdeModifier_3::*)(::out_stream &)) &ParabolicBoxDomainPdeModifier_3::OutputSimulationModifierParameters,
            " ", py::arg("rParamsFile"))
    ;
}
