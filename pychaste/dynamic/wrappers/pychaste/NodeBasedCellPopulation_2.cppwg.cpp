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
#include "AbstractCellBasedSimulation.hpp"
#include "AbstractImmersedBoundaryDivisionRule.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellDeltaNotchWriter.hpp"
#include "CellDivisionLocationsWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellRemovalLocationsWriter.hpp"
#include "CellRosetteRankWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "LegacyCellProliferativeTypesWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "PottsMeshWriter.hpp"
#include "RadialCellDataDistributionWriter.hpp"
#include "VertexIntersectionSwapLocationsWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "NodeBasedCellPopulation_2.cppwg.hpp"

namespace py = pybind11;
typedef NodeBasedCellPopulation<2> NodeBasedCellPopulation_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<2, 2> * _TetrahedralMesh_lt_2_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::vector<std::pair<Node<2> *, Node<2> *>> & _std_vector_lt_std_pair_lt_Node_lt_2_gt_Ptr_Node_lt_2_gt_Ptr_gt__gt_Ref;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;

class NodeBasedCellPopulation_2_Overrides : public NodeBasedCellPopulation_2
{
public:
    using NodeBasedCellPopulation_2::NodeBasedCellPopulation;
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> & rNewLocation) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            SetNode,
            nodeIndex,
            rNewLocation);
    }
    unsigned int GetNumNodes() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodeBasedCellPopulation_2,
            GetNumNodes,
            );
    }
    ::CellPtr GetCellUsingLocationIndex(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            NodeBasedCellPopulation_2,
            GetCellUsingLocationIndex,
            index);
    }
    ::Node<2> * GetNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            NodeBasedCellPopulation_2,
            GetNode,
            index);
    }
    unsigned int RemoveDeadCells() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            NodeBasedCellPopulation_2,
            RemoveDeadCells,
            );
    }
    void Update(bool hasHadBirthsOrDeaths) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            Update,
            hasHadBirthsOrDeaths);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>> pPopulationWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            AcceptPopulationWriter,
            pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>> pPopulationCountWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            AcceptPopulationCountWriter,
            pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>> pPopulationEventWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            AcceptPopulationEventWriter,
            pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<2, 2>> pCellWriter, ::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            AcceptCellWriter,
            pCellWriter,
            pCell);
    }
    double GetWidth(unsigned int const & rDimension) override
    {
        PYBIND11_OVERRIDE(
            double,
            NodeBasedCellPopulation_2,
            GetWidth,
            rDimension);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            NodeBasedCellPopulation_2,
            GetNeighbouringNodeIndices,
            index);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            NodeBasedCellPopulation_2,
            AddCell,
            pNewCell,
            pParentCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            NodeBasedCellPopulation_2,
            GetVolumeOfCell,
            pCell);
    }
    void UpdateCellProcessLocation() override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            UpdateCellProcessLocation,
            );
    }
    void UpdateParticlesAfterReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            UpdateParticlesAfterReMesh,
            rMap);
    }
    void Validate() override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulation_2,
            Validate,
            );
    }
};

void register_NodeBasedCellPopulation_2_class(py::module &m)
{
    py::class_<NodeBasedCellPopulation_2, NodeBasedCellPopulation_2_Overrides, boost::shared_ptr<NodeBasedCellPopulation_2>, AbstractCentreBasedCellPopulation<2>>(m, "NodeBasedCellPopulation_2")
        .def(py::init<::NodesOnlyMesh<2> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool, bool>(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("deleteMesh") = false, py::arg("validate") = true)
        .def(py::init<::NodesOnlyMesh<2> &>(), py::arg("rMesh"))
        .def("SetNode",
            (void(NodeBasedCellPopulation_2::*)(unsigned int, ::ChastePoint<2> &)) &NodeBasedCellPopulation_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("rNewLocation"))
        .def("GetNumNodes",
            (unsigned int(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::GetNumNodes,
            " ")
        .def("GetCellUsingLocationIndex",
            (::CellPtr(NodeBasedCellPopulation_2::*)(unsigned int)) &NodeBasedCellPopulation_2::GetCellUsingLocationIndex,
            " ", py::arg("index"))
        .def("GetNode",
            (::Node<2> *(NodeBasedCellPopulation_2::*)(unsigned int)) &NodeBasedCellPopulation_2::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("RemoveDeadCells",
            (unsigned int(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::RemoveDeadCells,
            " ")
        .def("Clear",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::Clear,
            " ")
        .def("Update",
            (void(NodeBasedCellPopulation_2::*)(bool)) &NodeBasedCellPopulation_2::Update,
            " ", py::arg("hasHadBirthsOrDeaths") = true)
        .def("OutputCellPopulationParameters",
            (void(NodeBasedCellPopulation_2::*)(::out_stream &)) &NodeBasedCellPopulation_2::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("AcceptPopulationWriter",
            (void(NodeBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>>)) &NodeBasedCellPopulation_2::AcceptPopulationWriter,
            " ", py::arg("pPopulationWriter"))
        .def("AcceptPopulationCountWriter",
            (void(NodeBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>>)) &NodeBasedCellPopulation_2::AcceptPopulationCountWriter,
            " ", py::arg("pPopulationCountWriter"))
        .def("AcceptPopulationEventWriter",
            (void(NodeBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>>)) &NodeBasedCellPopulation_2::AcceptPopulationEventWriter,
            " ", py::arg("pPopulationEventWriter"))
        .def("AcceptCellWriter",
            (void(NodeBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellWriter<2, 2>>, ::CellPtr)) &NodeBasedCellPopulation_2::AcceptCellWriter,
            " ", py::arg("pCellWriter"), py::arg("pCell"))
        .def("GetMechanicsCutOffLength",
            (double(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::GetMechanicsCutOffLength,
            " ")
        .def("GetUseVariableRadii",
            (bool(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::GetUseVariableRadii,
            " ")
        .def("SetUseVariableRadii",
            (void(NodeBasedCellPopulation_2::*)(bool)) &NodeBasedCellPopulation_2::SetUseVariableRadii,
            " ", py::arg("useVariableRadii") = true)
        .def("SetLoadBalanceMesh",
            (void(NodeBasedCellPopulation_2::*)(bool)) &NodeBasedCellPopulation_2::SetLoadBalanceMesh,
            " ", py::arg("loadBalanceMesh"))
        .def("SetLoadBalanceFrequency",
            (void(NodeBasedCellPopulation_2::*)(unsigned int)) &NodeBasedCellPopulation_2::SetLoadBalanceFrequency,
            " ", py::arg("loadBalanceFrequency"))
        .def("GetWidth",
            (double(NodeBasedCellPopulation_2::*)(unsigned int const &)) &NodeBasedCellPopulation_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("GetSizeOfCellPopulation",
            (::boost::numeric::ublas::c_vector<double, 2>(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::GetSizeOfCellPopulation,
            " ")
        .def("GetNodesWithinNeighbourhoodRadius",
            (::std::set<unsigned int>(NodeBasedCellPopulation_2::*)(unsigned int, double)) &NodeBasedCellPopulation_2::GetNodesWithinNeighbourhoodRadius,
            " ", py::arg("index"), py::arg("neighbourhoodRadius"))
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(NodeBasedCellPopulation_2::*)(unsigned int)) &NodeBasedCellPopulation_2::GetNeighbouringNodeIndices,
            " ", py::arg("index"))
        .def("AddCell",
            (::CellPtr(NodeBasedCellPopulation_2::*)(::CellPtr, ::CellPtr)) &NodeBasedCellPopulation_2::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell"))
        .def("GetVolumeOfCell",
            (double(NodeBasedCellPopulation_2::*)(::CellPtr)) &NodeBasedCellPopulation_2::GetVolumeOfCell,
            " ", py::arg("pCell"))
        .def("SendCellsToNeighbourProcesses",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::SendCellsToNeighbourProcesses,
            " ")
        .def("NonBlockingSendCellsToNeighbourProcesses",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::NonBlockingSendCellsToNeighbourProcesses,
            " ")
        .def("GetReceivedCells",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::GetReceivedCells,
            " ")
        .def("GetCellNodePair",
            (::std::pair<boost::shared_ptr<Cell>, Node<2> *>(NodeBasedCellPopulation_2::*)(unsigned int)) &NodeBasedCellPopulation_2::GetCellNodePair,
            " ", py::arg("nodeIndex"))
        .def("AddReceivedCells",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::AddReceivedCells,
            " ")
        .def("UpdateCellProcessLocation",
            (void(NodeBasedCellPopulation_2::*)()) &NodeBasedCellPopulation_2::UpdateCellProcessLocation,
            " ")
        .def("AddPopulationWriterVoronoiDataWriter", &NodeBasedCellPopulation_2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &NodeBasedCellPopulation_2::AddCellWriter<CellLabelWriter>)
    ;
}
