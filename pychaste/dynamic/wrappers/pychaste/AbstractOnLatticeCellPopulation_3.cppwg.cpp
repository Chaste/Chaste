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
#include "AbstractOnLatticeCellPopulation.hpp"

#include "AbstractOnLatticeCellPopulation_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOnLatticeCellPopulation<3> AbstractOnLatticeCellPopulation_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const;

class AbstractOnLatticeCellPopulation_3_Overrides : public AbstractOnLatticeCellPopulation_3
{
public:
    using AbstractOnLatticeCellPopulation_3::AbstractOnLatticeCellPopulation;
    void UpdateCellLocations(double dt) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation_3,
            UpdateCellLocations,
            dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<3> & rNewLocation) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation_3,
            SetNode,
            index,
            rNewLocation);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            AbstractOnLatticeCellPopulation_3,
            GetNeighbouringNodeIndices,
            index);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation_3,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    double GetDefaultTimeStep() override
    {
        PYBIND11_OVERRIDE(
            double,
            AbstractOnLatticeCellPopulation_3,
            GetDefaultTimeStep,
            );
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<3>> pUpdateRule) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation_3,
            AddUpdateRule,
            pUpdateRule);
    }
    void RemoveAllUpdateRules() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation_3,
            RemoveAllUpdateRules,
            );
    }
    ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const GetUpdateRuleCollection() const override
    {
        PYBIND11_OVERRIDE(
            _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const,
            AbstractOnLatticeCellPopulation_3,
            GetUpdateRuleCollection,
            );
    }
};

void register_AbstractOnLatticeCellPopulation_3_class(py::module &m)
{
    py::class_<AbstractOnLatticeCellPopulation_3, AbstractOnLatticeCellPopulation_3_Overrides, boost::shared_ptr<AbstractOnLatticeCellPopulation_3>, AbstractCellPopulation<3>>(m, "AbstractOnLatticeCellPopulation_3")
        .def("UpdateCellLocations",
            (void(AbstractOnLatticeCellPopulation_3::*)(double)) &AbstractOnLatticeCellPopulation_3::UpdateCellLocations,
            " ", py::arg("dt"))
        .def("GetUpdateNodesInRandomOrder",
            (bool(AbstractOnLatticeCellPopulation_3::*)()) &AbstractOnLatticeCellPopulation_3::GetUpdateNodesInRandomOrder,
            " ")
        .def("SetUpdateNodesInRandomOrder",
            (void(AbstractOnLatticeCellPopulation_3::*)(bool)) &AbstractOnLatticeCellPopulation_3::SetUpdateNodesInRandomOrder,
            " ", py::arg("updateNodesInRandomOrder"))
        .def("SetIterateRandomlyOverUpdateRuleCollection",
            (void(AbstractOnLatticeCellPopulation_3::*)(bool)) &AbstractOnLatticeCellPopulation_3::SetIterateRandomlyOverUpdateRuleCollection,
            " ", py::arg("iterateRandomly"))
        .def("GetIterateRandomlyOverUpdateRuleCollection",
            (bool(AbstractOnLatticeCellPopulation_3::*)()) &AbstractOnLatticeCellPopulation_3::GetIterateRandomlyOverUpdateRuleCollection,
            " ")
        .def("SetNode",
            (void(AbstractOnLatticeCellPopulation_3::*)(unsigned int, ::ChastePoint<3> &)) &AbstractOnLatticeCellPopulation_3::SetNode,
            " ", py::arg("index"), py::arg("rNewLocation"))
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(AbstractOnLatticeCellPopulation_3::*)(unsigned int)) &AbstractOnLatticeCellPopulation_3::GetNeighbouringNodeIndices,
            " ", py::arg("index"))
        .def("OutputCellPopulationParameters",
            (void(AbstractOnLatticeCellPopulation_3::*)(::out_stream &)) &AbstractOnLatticeCellPopulation_3::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("GetDefaultTimeStep",
            (double(AbstractOnLatticeCellPopulation_3::*)()) &AbstractOnLatticeCellPopulation_3::GetDefaultTimeStep,
            " ")
        .def("AddUpdateRule",
            (void(AbstractOnLatticeCellPopulation_3::*)(::boost::shared_ptr<AbstractUpdateRule<3>>)) &AbstractOnLatticeCellPopulation_3::AddUpdateRule,
            " ", py::arg("pUpdateRule"))
        .def("RemoveAllUpdateRules",
            (void(AbstractOnLatticeCellPopulation_3::*)()) &AbstractOnLatticeCellPopulation_3::RemoveAllUpdateRules,
            " ")
        .def("GetUpdateRuleCollection",
            (::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const(AbstractOnLatticeCellPopulation_3::*)() const) &AbstractOnLatticeCellPopulation_3::GetUpdateRuleCollection,
            " ")
    ;
}
