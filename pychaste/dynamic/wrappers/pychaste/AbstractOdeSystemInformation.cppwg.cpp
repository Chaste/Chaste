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
#include "AbstractOdeSystemInformation.hpp"

#include "AbstractOdeSystemInformation.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeSystemInformation AbstractOdeSystemInformation;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeSystemInformation_Overrides : public AbstractOdeSystemInformation
{
public:
    using AbstractOdeSystemInformation::AbstractOdeSystemInformation;
    void Initialise() override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOdeSystemInformation,
            Initialise,
            );
    }
};

void register_AbstractOdeSystemInformation_class(py::module &m)
{
    py::class_<AbstractOdeSystemInformation, AbstractOdeSystemInformation_Overrides, boost::shared_ptr<AbstractOdeSystemInformation>>(m, "AbstractOdeSystemInformation")
        .def(py::init<>())
        .def("GetSystemName",
            (::std::string(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetSystemName,
            " ")
        .def("GetFreeVariableName",
            (::std::string(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetFreeVariableName,
            " ")
        .def("GetFreeVariableUnits",
            (::std::string(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetFreeVariableUnits,
            " ")
        .def("SetDefaultInitialConditions",
            (void(AbstractOdeSystemInformation::*)(::std::vector<double> const &)) &AbstractOdeSystemInformation::SetDefaultInitialConditions,
            " ", py::arg("rInitialConditions"))
        .def("SetDefaultInitialCondition",
            (void(AbstractOdeSystemInformation::*)(unsigned int, double)) &AbstractOdeSystemInformation::SetDefaultInitialCondition,
            " ", py::arg("index"), py::arg("initialCondition"))
        .def("GetInitialConditions",
            (::std::vector<double>(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetInitialConditions,
            " ")
        .def("rGetStateVariableNames",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetStateVariableNames,
            " ", py::return_value_policy::reference_internal)
        .def("rGetStateVariableUnits",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetStateVariableUnits,
            " ", py::return_value_policy::reference_internal)
        .def("GetStateVariableIndex",
            (unsigned int(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::GetStateVariableIndex,
            " ", py::arg("rName"))
        .def("HasStateVariable",
            (bool(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::HasStateVariable,
            " ", py::arg("rName"))
        .def("GetStateVariableUnits",
            (::std::string(AbstractOdeSystemInformation::*)(unsigned int) const) &AbstractOdeSystemInformation::GetStateVariableUnits,
            " ", py::arg("index"))
        .def("rGetParameterNames",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetParameterNames,
            " ", py::return_value_policy::reference_internal)
        .def("rGetParameterUnits",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetParameterUnits,
            " ", py::return_value_policy::reference_internal)
        .def("GetParameterIndex",
            (unsigned int(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::GetParameterIndex,
            " ", py::arg("rName"))
        .def("HasParameter",
            (bool(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::HasParameter,
            " ", py::arg("rName"))
        .def("GetParameterUnits",
            (::std::string(AbstractOdeSystemInformation::*)(unsigned int) const) &AbstractOdeSystemInformation::GetParameterUnits,
            " ", py::arg("index"))
        .def("GetNumberOfParameters",
            (unsigned int(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetNumberOfParameters,
            " ")
        .def("GetAnyVariableIndex",
            (unsigned int(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::GetAnyVariableIndex,
            " ", py::arg("rName"))
        .def("HasAnyVariable",
            (bool(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::HasAnyVariable,
            " ", py::arg("rName"))
        .def("GetAnyVariableUnits",
            (::std::string(AbstractOdeSystemInformation::*)(unsigned int) const) &AbstractOdeSystemInformation::GetAnyVariableUnits,
            " ", py::arg("index"))
        .def("rGetDerivedQuantityNames",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetDerivedQuantityNames,
            " ", py::return_value_policy::reference_internal)
        .def("rGetDerivedQuantityUnits",
            (::std::vector<std::basic_string<char>> const &(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::rGetDerivedQuantityUnits,
            " ", py::return_value_policy::reference_internal)
        .def("GetDerivedQuantityIndex",
            (unsigned int(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::GetDerivedQuantityIndex,
            " ", py::arg("rName"))
        .def("HasDerivedQuantity",
            (bool(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::HasDerivedQuantity,
            " ", py::arg("rName"))
        .def("GetDerivedQuantityUnits",
            (::std::string(AbstractOdeSystemInformation::*)(unsigned int) const) &AbstractOdeSystemInformation::GetDerivedQuantityUnits,
            " ", py::arg("index"))
        .def("GetNumberOfDerivedQuantities",
            (unsigned int(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetNumberOfDerivedQuantities,
            " ")
        .def("GetNumberOfAttributes",
            (unsigned int(AbstractOdeSystemInformation::*)() const) &AbstractOdeSystemInformation::GetNumberOfAttributes,
            " ")
        .def("HasAttribute",
            (bool(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::HasAttribute,
            " ", py::arg("rName"))
        .def("GetAttribute",
            (double(AbstractOdeSystemInformation::*)(::std::string const &) const) &AbstractOdeSystemInformation::GetAttribute,
            " ", py::arg("rName"))
    ;
}
