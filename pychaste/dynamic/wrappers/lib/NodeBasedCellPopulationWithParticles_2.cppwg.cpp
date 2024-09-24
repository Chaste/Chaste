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
#include "NodeBasedCellPopulationWithParticles.hpp"

#include "NodeBasedCellPopulationWithParticles_2.cppwg.hpp"

namespace py = pybind11;
typedef NodeBasedCellPopulationWithParticles<2> NodeBasedCellPopulationWithParticles_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::CellPtr _CellPtr;

class NodeBasedCellPopulationWithParticles_2_Overrides : public NodeBasedCellPopulationWithParticles_2
{
public:
    using NodeBasedCellPopulationWithParticles_2::NodeBasedCellPopulationWithParticles;
    void UpdateParticlesAfterReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles_2,
            UpdateParticlesAfterReMesh,
            rMap);
    }
    bool IsParticle(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            bool,
            NodeBasedCellPopulationWithParticles_2,
            IsParticle,
            index);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            NodeBasedCellPopulationWithParticles_2,
            AddCell,
            pNewCell,
            pParentCell);
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles_2,
            WriteVtkResultsToFile,
            rDirectory);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles_2,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    void Validate() override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles_2,
            Validate,
            );
    }
    void AcceptCellWritersAcrossPopulation() override
    {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles_2,
            AcceptCellWritersAcrossPopulation,
            );
    }
};

void register_NodeBasedCellPopulationWithParticles_2_class(py::module &m)
{
    py::class_<NodeBasedCellPopulationWithParticles_2, NodeBasedCellPopulationWithParticles_2_Overrides, boost::shared_ptr<NodeBasedCellPopulationWithParticles_2>, NodeBasedCellPopulation<2>>(m, "NodeBasedCellPopulationWithParticles_2")
        .def(py::init<::NodesOnlyMesh<2> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool>(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("deleteMesh") = false)
        .def(py::init<::NodesOnlyMesh<2> &>(), py::arg("rMesh"))
        .def("UpdateParticlesAfterReMesh",
            (void(NodeBasedCellPopulationWithParticles_2::*)(::NodeMap &)) &NodeBasedCellPopulationWithParticles_2::UpdateParticlesAfterReMesh,
            " ", py::arg("rMap"))
        .def("IsParticle",
            (bool(NodeBasedCellPopulationWithParticles_2::*)(unsigned int)) &NodeBasedCellPopulationWithParticles_2::IsParticle,
            " ", py::arg("index"))
        .def("GetParticleIndices",
            (::std::set<unsigned int>(NodeBasedCellPopulationWithParticles_2::*)()) &NodeBasedCellPopulationWithParticles_2::GetParticleIndices,
            " ")
        .def("AddCell",
            (::CellPtr(NodeBasedCellPopulationWithParticles_2::*)(::CellPtr, ::CellPtr)) &NodeBasedCellPopulationWithParticles_2::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell"))
        .def("WriteVtkResultsToFile",
            (void(NodeBasedCellPopulationWithParticles_2::*)(::std::string const &)) &NodeBasedCellPopulationWithParticles_2::WriteVtkResultsToFile,
            " ", py::arg("rDirectory"))
        .def("OutputCellPopulationParameters",
            (void(NodeBasedCellPopulationWithParticles_2::*)(::out_stream &)) &NodeBasedCellPopulationWithParticles_2::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("AddPopulationWriterVoronoiDataWriter", &NodeBasedCellPopulationWithParticles_2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &NodeBasedCellPopulationWithParticles_2::AddCellWriter<CellLabelWriter>)
    ;
}
