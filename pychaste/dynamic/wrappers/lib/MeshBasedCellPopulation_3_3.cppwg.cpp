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
#include "MeshBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation_3_3.cppwg.hpp"

namespace py = pybind11;
typedef MeshBasedCellPopulation<3, 3> MeshBasedCellPopulation_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::vector<std::pair<Node<3> *, Node<3> *>> & _std_vector_lt_std_pair_lt_Node_lt_3_gt_Ptr_Node_lt_3_gt_Ptr_gt__gt_Ref;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;

class MeshBasedCellPopulation_3_3_Overrides : public MeshBasedCellPopulation_3_3
{
public:
    using MeshBasedCellPopulation_3_3::MeshBasedCellPopulation;
    unsigned int AddNode(::Node<3> * pNewNode) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_3_3,
            AddNode,
            pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> & rNewLocation) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            SetNode,
            nodeIndex,
            rNewLocation);
    }
    double GetDampingConstant(unsigned int nodeIndex) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_3_3,
            GetDampingConstant,
            nodeIndex);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            OpenWritersFiles,
            rOutputFileHandler);
    }
    unsigned int RemoveDeadCells() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_3_3,
            RemoveDeadCells,
            );
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            MeshBasedCellPopulation_3_3,
            AddCell,
            pNewCell,
            pParentCell);
    }
    void WriteResultsToFiles(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            WriteResultsToFiles,
            rDirectory);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>> pPopulationWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            AcceptPopulationWriter,
            pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>> pPopulationCountWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            AcceptPopulationCountWriter,
            pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>> pPopulationEventWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            AcceptPopulationEventWriter,
            pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<3, 3>> pCellWriter, ::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            AcceptCellWriter,
            pCellWriter,
            pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            Update,
            hasHadBirthsOrDeaths);
    }
    ::Node<3> * GetNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            MeshBasedCellPopulation_3_3,
            GetNode,
            index);
    }
    unsigned int GetNumNodes() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_3_3,
            GetNumNodes,
            );
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            WriteVtkResultsToFile,
            rDirectory);
    }
    double GetVolumeOfCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_3_3,
            GetVolumeOfCell,
            pCell);
    }
    double GetWidth(unsigned int const & rDimension) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_3_3,
            GetWidth,
            rDimension);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            WriteDataToVisualizerSetupFile,
            pVizSetupFile);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            MeshBasedCellPopulation_3_3,
            GetNeighbouringNodeIndices,
            index);
    }
    void UpdateGhostNodesAfterReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            UpdateGhostNodesAfterReMesh,
            rMap);
    }
    void Validate() override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_3_3,
            Validate,
            );
    }
};

void register_MeshBasedCellPopulation_3_3_class(py::module &m)
{
    py::class_<MeshBasedCellPopulation_3_3, MeshBasedCellPopulation_3_3_Overrides, boost::shared_ptr<MeshBasedCellPopulation_3_3>, AbstractCentreBasedCellPopulation<3>>(m, "MeshBasedCellPopulation_3_3")
        .def(py::init<::MutableMesh<3, 3> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool, bool>(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = ::std::vector<unsigned int> {}, py::arg("deleteMesh") = false, py::arg("validate") = true)
        .def(py::init<::MutableMesh<3, 3> &>(), py::arg("rMesh"))
        .def("UseAreaBasedDampingConstant",
            (bool(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::UseAreaBasedDampingConstant,
            " ")
        .def("AddNode",
            (unsigned int(MeshBasedCellPopulation_3_3::*)(::Node<3> *)) &MeshBasedCellPopulation_3_3::AddNode,
            " ", py::arg("pNewNode"))
        .def("SetNode",
            (void(MeshBasedCellPopulation_3_3::*)(unsigned int, ::ChastePoint<3> &)) &MeshBasedCellPopulation_3_3::SetNode,
            " ", py::arg("nodeIndex"), py::arg("rNewLocation"))
        .def("GetDampingConstant",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int)) &MeshBasedCellPopulation_3_3::GetDampingConstant,
            " ", py::arg("nodeIndex"))
        .def("SetAreaBasedDampingConstant",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::SetAreaBasedDampingConstant,
            " ", py::arg("useAreaBasedDampingConstant"))
        .def("OpenWritersFiles",
            (void(MeshBasedCellPopulation_3_3::*)(::OutputFileHandler &)) &MeshBasedCellPopulation_3_3::OpenWritersFiles,
            " ", py::arg("rOutputFileHandler"))
        .def("RemoveDeadCells",
            (unsigned int(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::RemoveDeadCells,
            " ")
        .def("AddCell",
            (::CellPtr(MeshBasedCellPopulation_3_3::*)(::CellPtr, ::CellPtr)) &MeshBasedCellPopulation_3_3::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell"))
        .def("WriteResultsToFiles",
            (void(MeshBasedCellPopulation_3_3::*)(::std::string const &)) &MeshBasedCellPopulation_3_3::WriteResultsToFiles,
            " ", py::arg("rDirectory"))
        .def("AcceptPopulationWriter",
            (void(MeshBasedCellPopulation_3_3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &MeshBasedCellPopulation_3_3::AcceptPopulationWriter,
            " ", py::arg("pPopulationWriter"))
        .def("AcceptPopulationCountWriter",
            (void(MeshBasedCellPopulation_3_3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &MeshBasedCellPopulation_3_3::AcceptPopulationCountWriter,
            " ", py::arg("pPopulationCountWriter"))
        .def("AcceptPopulationEventWriter",
            (void(MeshBasedCellPopulation_3_3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &MeshBasedCellPopulation_3_3::AcceptPopulationEventWriter,
            " ", py::arg("pPopulationEventWriter"))
        .def("AcceptCellWriter",
            (void(MeshBasedCellPopulation_3_3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>, ::CellPtr)) &MeshBasedCellPopulation_3_3::AcceptCellWriter,
            " ", py::arg("pCellWriter"), py::arg("pCell"))
        .def("Update",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::Update,
            " ", py::arg("hasHadBirthsOrDeaths") = true)
        .def("TessellateIfNeeded",
            (void(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::TessellateIfNeeded,
            " ")
        .def("DivideLongSprings",
            (void(MeshBasedCellPopulation_3_3::*)(double)) &MeshBasedCellPopulation_3_3::DivideLongSprings,
            " ", py::arg("springDivisionThreshold"))
        .def("GetNode",
            (::Node<3> *(MeshBasedCellPopulation_3_3::*)(unsigned int)) &MeshBasedCellPopulation_3_3::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNumNodes",
            (unsigned int(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetNumNodes,
            " ")
        .def("WriteVtkResultsToFile",
            (void(MeshBasedCellPopulation_3_3::*)(::std::string const &)) &MeshBasedCellPopulation_3_3::WriteVtkResultsToFile,
            " ", py::arg("rDirectory"))
        .def("GetVolumeOfCell",
            (double(MeshBasedCellPopulation_3_3::*)(::CellPtr)) &MeshBasedCellPopulation_3_3::GetVolumeOfCell,
            " ", py::arg("pCell"))
        .def("CreateVoronoiTessellation",
            (void(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::CreateVoronoiTessellation,
            " ")
        .def("GetVolumeOfVoronoiElement",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int)) &MeshBasedCellPopulation_3_3::GetVolumeOfVoronoiElement,
            " ", py::arg("index"))
        .def("GetSurfaceAreaOfVoronoiElement",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int)) &MeshBasedCellPopulation_3_3::GetSurfaceAreaOfVoronoiElement,
            " ", py::arg("index"))
        .def("GetVoronoiEdgeLength",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int, unsigned int)) &MeshBasedCellPopulation_3_3::GetVoronoiEdgeLength,
            " ", py::arg("index1"), py::arg("index2"))
        .def("GetWidth",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int const &)) &MeshBasedCellPopulation_3_3::GetWidth,
            " ", py::arg("rDimension"))
        .def("WriteDataToVisualizerSetupFile",
            (void(MeshBasedCellPopulation_3_3::*)(::out_stream &)) &MeshBasedCellPopulation_3_3::WriteDataToVisualizerSetupFile,
            " ", py::arg("pVizSetupFile"))
        .def("SpringsBegin",
            (::MeshBasedCellPopulation<3>::SpringIterator(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::SpringsBegin,
            " ")
        .def("SpringsEnd",
            (::MeshBasedCellPopulation<3>::SpringIterator(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::SpringsEnd,
            " ")
        .def("CheckCellPointers",
            (void(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::CheckCellPointers,
            " ")
        .def("GetAreaBasedDampingConstantParameter",
            (double(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetAreaBasedDampingConstantParameter,
            " ")
        .def("SetAreaBasedDampingConstantParameter",
            (void(MeshBasedCellPopulation_3_3::*)(double)) &MeshBasedCellPopulation_3_3::SetAreaBasedDampingConstantParameter,
            " ", py::arg("areaBasedDampingConstantParameter"))
        .def("OutputCellPopulationParameters",
            (void(MeshBasedCellPopulation_3_3::*)(::out_stream &)) &MeshBasedCellPopulation_3_3::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("SetWriteVtkAsPoints",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::SetWriteVtkAsPoints,
            " ", py::arg("writeVtkAsPoints"))
        .def("GetWriteVtkAsPoints",
            (bool(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetWriteVtkAsPoints,
            " ")
        .def("SetBoundVoronoiTessellation",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::SetBoundVoronoiTessellation,
            " ", py::arg("boundVoronoiTessellation"))
        .def("GetBoundVoronoiTessellation",
            (bool(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetBoundVoronoiTessellation,
            " ")
        .def("SetScaleBoundByEdgeLength",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::SetScaleBoundByEdgeLength,
            " ", py::arg("scaleBoundByEdgeLength"))
        .def("GetScaleBoundByEdgeLength",
            (bool(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetScaleBoundByEdgeLength,
            " ")
        .def("SetBoundedVoroniTesselationLengthCutoff",
            (void(MeshBasedCellPopulation_3_3::*)(double)) &MeshBasedCellPopulation_3_3::SetBoundedVoroniTesselationLengthCutoff,
            " ", py::arg("boundedVoroniTesselationLengthCutoff"))
        .def("GetBoundedVoroniTesselationLengthCutoff",
            (double(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetBoundedVoroniTesselationLengthCutoff,
            " ")
        .def("SetOffsetNewBoundaryNodes",
            (void(MeshBasedCellPopulation_3_3::*)(bool)) &MeshBasedCellPopulation_3_3::SetOffsetNewBoundaryNodes,
            " ", py::arg("offsetNewBoundaryNodes"))
        .def("GetOffsetNewBoundaryNodes",
            (bool(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::GetOffsetNewBoundaryNodes,
            " ")
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(MeshBasedCellPopulation_3_3::*)(unsigned int)) &MeshBasedCellPopulation_3_3::GetNeighbouringNodeIndices,
            " ", py::arg("index"))
        .def("CalculateRestLengths",
            (void(MeshBasedCellPopulation_3_3::*)()) &MeshBasedCellPopulation_3_3::CalculateRestLengths,
            " ")
        .def("GetRestLength",
            (double(MeshBasedCellPopulation_3_3::*)(unsigned int, unsigned int)) &MeshBasedCellPopulation_3_3::GetRestLength,
            " ", py::arg("indexA"), py::arg("indexB"))
        .def("SetRestLength",
            (void(MeshBasedCellPopulation_3_3::*)(unsigned int, unsigned int, double)) &MeshBasedCellPopulation_3_3::SetRestLength,
            " ", py::arg("indexA"), py::arg("indexB"), py::arg("restLength"))
        .def("AddPopulationWriterVoronoiDataWriter", &MeshBasedCellPopulation_3_3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &MeshBasedCellPopulation_3_3::AddCellWriter<CellLabelWriter>)
    ;
}
