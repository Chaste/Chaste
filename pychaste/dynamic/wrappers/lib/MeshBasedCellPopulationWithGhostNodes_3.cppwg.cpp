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
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "MeshBasedCellPopulationWithGhostNodes_3.cppwg.hpp"

namespace py = pybind11;
typedef MeshBasedCellPopulationWithGhostNodes<3> MeshBasedCellPopulationWithGhostNodes_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::CellPtr _CellPtr;

class MeshBasedCellPopulationWithGhostNodes_3_Overrides : public MeshBasedCellPopulationWithGhostNodes_3
{
public:
    using MeshBasedCellPopulationWithGhostNodes_3::MeshBasedCellPopulationWithGhostNodes;
    ::TetrahedralMesh<3, 3> * GetTetrahedralMeshForPdeModifier() override
    {
        PYBIND11_OVERRIDE(
            _TetrahedralMesh_lt_3_3_gt_Ptr,
            MeshBasedCellPopulationWithGhostNodes_3,
            GetTetrahedralMeshForPdeModifier,
            );
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            MeshBasedCellPopulationWithGhostNodes_3,
            GetNeighbouringLocationIndices,
            pCell);
    }
    bool IsGhostNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            bool,
            MeshBasedCellPopulationWithGhostNodes_3,
            IsGhostNode,
            index);
    }
    void UpdateGhostNodesAfterReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulationWithGhostNodes_3,
            UpdateGhostNodesAfterReMesh,
            rMap);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            MeshBasedCellPopulationWithGhostNodes_3,
            AddCell,
            pNewCell,
            pParentCell);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulationWithGhostNodes_3,
            OpenWritersFiles,
            rOutputFileHandler);
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulationWithGhostNodes_3,
            WriteVtkResultsToFile,
            rDirectory);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulationWithGhostNodes_3,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    void AcceptCellWritersAcrossPopulation() override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulationWithGhostNodes_3,
            AcceptCellWritersAcrossPopulation,
            );
    }
};

void register_MeshBasedCellPopulationWithGhostNodes_3_class(py::module &m)
{
    py::class_<MeshBasedCellPopulationWithGhostNodes_3, MeshBasedCellPopulationWithGhostNodes_3_Overrides, boost::shared_ptr<MeshBasedCellPopulationWithGhostNodes_3>, MeshBasedCellPopulation<3>>(m, "MeshBasedCellPopulationWithGhostNodes_3")
        .def(py::init<::MutableMesh<3, 3> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool, double, double, double>(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("deleteMesh") = false, py::arg("ghostCellSpringStiffness") = 15., py::arg("ghostGhostSpringStiffness") = 15., py::arg("ghostSpringRestLength") = 1.)
        .def(py::init<::MutableMesh<3, 3> &, double, double, double>(), py::arg("rMesh"), py::arg("ghostCellSpringStiffness") = 15., py::arg("ghostGhostSpringStiffness") = 15., py::arg("ghostSpringRestLength") = 1.)
        .def("GetTetrahedralMeshForPdeModifier",
            (::TetrahedralMesh<3, 3> *(MeshBasedCellPopulationWithGhostNodes_3::*)()) &MeshBasedCellPopulationWithGhostNodes_3::GetTetrahedralMeshForPdeModifier,
            " ", py::return_value_policy::reference)
        .def("GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(MeshBasedCellPopulationWithGhostNodes_3::*)(::CellPtr)) &MeshBasedCellPopulationWithGhostNodes_3::GetNeighbouringLocationIndices,
            " ", py::arg("pCell"))
        .def("SetGhostNodes",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(::std::set<unsigned int> const &)) &MeshBasedCellPopulationWithGhostNodes_3::SetGhostNodes,
            " ", py::arg("rGhostNodeIndices"))
        .def("ApplyGhostForces",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)()) &MeshBasedCellPopulationWithGhostNodes_3::ApplyGhostForces,
            " ")
        .def("rGetGhostNodes",
            (::std::vector<bool> &(MeshBasedCellPopulationWithGhostNodes_3::*)()) &MeshBasedCellPopulationWithGhostNodes_3::rGetGhostNodes,
            " ", py::return_value_policy::reference_internal)
        .def("IsGhostNode",
            (bool(MeshBasedCellPopulationWithGhostNodes_3::*)(unsigned int)) &MeshBasedCellPopulationWithGhostNodes_3::IsGhostNode,
            " ", py::arg("index"))
        .def("GetGhostNodeIndices",
            (::std::set<unsigned int>(MeshBasedCellPopulationWithGhostNodes_3::*)()) &MeshBasedCellPopulationWithGhostNodes_3::GetGhostNodeIndices,
            " ")
        .def("RemoveGhostNode",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(unsigned int)) &MeshBasedCellPopulationWithGhostNodes_3::RemoveGhostNode,
            " ", py::arg("nodeIndex"))
        .def("UpdateGhostNodesAfterReMesh",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(::NodeMap &)) &MeshBasedCellPopulationWithGhostNodes_3::UpdateGhostNodesAfterReMesh,
            " ", py::arg("rMap"))
        .def("CalculateForceBetweenGhostNodes",
            (::boost::numeric::ublas::c_vector<double, 3>(MeshBasedCellPopulationWithGhostNodes_3::*)(unsigned int const &, unsigned int const &)) &MeshBasedCellPopulationWithGhostNodes_3::CalculateForceBetweenGhostNodes,
            " ", py::arg("rNodeAGlobalIndex"), py::arg("rNodeBGlobalIndex"))
        .def("AddCell",
            (::CellPtr(MeshBasedCellPopulationWithGhostNodes_3::*)(::CellPtr, ::CellPtr)) &MeshBasedCellPopulationWithGhostNodes_3::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell"))
        .def("OpenWritersFiles",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(::OutputFileHandler &)) &MeshBasedCellPopulationWithGhostNodes_3::OpenWritersFiles,
            " ", py::arg("rOutputFileHandler"))
        .def("WriteVtkResultsToFile",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(::std::string const &)) &MeshBasedCellPopulationWithGhostNodes_3::WriteVtkResultsToFile,
            " ", py::arg("rDirectory"))
        .def("OutputCellPopulationParameters",
            (void(MeshBasedCellPopulationWithGhostNodes_3::*)(::out_stream &)) &MeshBasedCellPopulationWithGhostNodes_3::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("AddPopulationWriterVoronoiDataWriter", &MeshBasedCellPopulationWithGhostNodes_3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &MeshBasedCellPopulationWithGhostNodes_3::AddCellWriter<CellLabelWriter>)
    ;
}
