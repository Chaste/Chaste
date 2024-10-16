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
#include "PythonUblasObjectConverters.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation_2_2.cppwg.hpp"

namespace py = pybind11;
typedef MeshBasedCellPopulation<2, 2> MeshBasedCellPopulation_2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<2, 2> * _TetrahedralMesh_lt_2_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::vector<std::pair<Node<2> *, Node<2> *>> & _std_vector_lt_std_pair_lt_Node_lt_2_gt_Ptr_Node_lt_2_gt_Ptr_gt__gt_Ref;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;

class MeshBasedCellPopulation_2_2_Overrides : public MeshBasedCellPopulation_2_2
{
public:
    using MeshBasedCellPopulation_2_2::MeshBasedCellPopulation;
    unsigned int AddNode(::Node<2> * pNewNode) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_2_2,
            AddNode,
            pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> & rNewLocation) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            SetNode,
            nodeIndex,
            rNewLocation);
    }
    double GetDampingConstant(unsigned int nodeIndex) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_2_2,
            GetDampingConstant,
            nodeIndex);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            OpenWritersFiles,
            rOutputFileHandler);
    }
    unsigned int RemoveDeadCells() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_2_2,
            RemoveDeadCells,
            );
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            MeshBasedCellPopulation_2_2,
            AddCell,
            pNewCell,
            pParentCell);
    }
    void WriteResultsToFiles(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            WriteResultsToFiles,
            rDirectory);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>> pPopulationWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            AcceptPopulationWriter,
            pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>> pPopulationCountWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            AcceptPopulationCountWriter,
            pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>> pPopulationEventWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            AcceptPopulationEventWriter,
            pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<2, 2>> pCellWriter, ::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            AcceptCellWriter,
            pCellWriter,
            pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            Update,
            hasHadBirthsOrDeaths);
    }
    ::Node<2> * GetNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            MeshBasedCellPopulation_2_2,
            GetNode,
            index);
    }
    unsigned int GetNumNodes() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            MeshBasedCellPopulation_2_2,
            GetNumNodes,
            );
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            WriteVtkResultsToFile,
            rDirectory);
    }
    double GetVolumeOfCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_2_2,
            GetVolumeOfCell,
            pCell);
    }
    double GetWidth(unsigned int const & rDimension) override
    {
        PYBIND11_OVERRIDE(
            double,
            MeshBasedCellPopulation_2_2,
            GetWidth,
            rDimension);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            WriteDataToVisualizerSetupFile,
            pVizSetupFile);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            MeshBasedCellPopulation_2_2,
            GetNeighbouringNodeIndices,
            index);
    }
    void UpdateGhostNodesAfterReMesh(::NodeMap & rMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            UpdateGhostNodesAfterReMesh,
            rMap);
    }
    void Validate() override
    {
        PYBIND11_OVERRIDE(
            void,
            MeshBasedCellPopulation_2_2,
            Validate,
            );
    }
};

void register_MeshBasedCellPopulation_2_2_class(py::module &m)
{
    py::class_<MeshBasedCellPopulation_2_2, MeshBasedCellPopulation_2_2_Overrides, boost::shared_ptr<MeshBasedCellPopulation_2_2>, AbstractCentreBasedCellPopulation<2>>(m, "MeshBasedCellPopulation_2_2")
        .def(py::init<::MutableMesh<2, 2> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool, bool>(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = ::std::vector<unsigned int> {}, py::arg("deleteMesh") = false, py::arg("validate") = true)
        .def(py::init<::MutableMesh<2, 2> &>(), py::arg("rMesh"))
        .def("UseAreaBasedDampingConstant",
            (bool(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::UseAreaBasedDampingConstant,
            " ")
        .def("AddNode",
            (unsigned int(MeshBasedCellPopulation_2_2::*)(::Node<2> *)) &MeshBasedCellPopulation_2_2::AddNode,
            " ", py::arg("pNewNode"))
        .def("SetNode",
            (void(MeshBasedCellPopulation_2_2::*)(unsigned int, ::ChastePoint<2> &)) &MeshBasedCellPopulation_2_2::SetNode,
            " ", py::arg("nodeIndex"), py::arg("rNewLocation"))
        .def("GetDampingConstant",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int)) &MeshBasedCellPopulation_2_2::GetDampingConstant,
            " ", py::arg("nodeIndex"))
        .def("SetAreaBasedDampingConstant",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::SetAreaBasedDampingConstant,
            " ", py::arg("useAreaBasedDampingConstant"))
        .def("OpenWritersFiles",
            (void(MeshBasedCellPopulation_2_2::*)(::OutputFileHandler &)) &MeshBasedCellPopulation_2_2::OpenWritersFiles,
            " ", py::arg("rOutputFileHandler"))
        .def("RemoveDeadCells",
            (unsigned int(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::RemoveDeadCells,
            " ")
        .def("AddCell",
            (::CellPtr(MeshBasedCellPopulation_2_2::*)(::CellPtr, ::CellPtr)) &MeshBasedCellPopulation_2_2::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell"))
        .def("WriteResultsToFiles",
            (void(MeshBasedCellPopulation_2_2::*)(::std::string const &)) &MeshBasedCellPopulation_2_2::WriteResultsToFiles,
            " ", py::arg("rDirectory"))
        .def("AcceptPopulationWriter",
            (void(MeshBasedCellPopulation_2_2::*)(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>>)) &MeshBasedCellPopulation_2_2::AcceptPopulationWriter,
            " ", py::arg("pPopulationWriter"))
        .def("AcceptPopulationCountWriter",
            (void(MeshBasedCellPopulation_2_2::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>>)) &MeshBasedCellPopulation_2_2::AcceptPopulationCountWriter,
            " ", py::arg("pPopulationCountWriter"))
        .def("AcceptPopulationEventWriter",
            (void(MeshBasedCellPopulation_2_2::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>>)) &MeshBasedCellPopulation_2_2::AcceptPopulationEventWriter,
            " ", py::arg("pPopulationEventWriter"))
        .def("AcceptCellWriter",
            (void(MeshBasedCellPopulation_2_2::*)(::boost::shared_ptr<AbstractCellWriter<2, 2>>, ::CellPtr)) &MeshBasedCellPopulation_2_2::AcceptCellWriter,
            " ", py::arg("pCellWriter"), py::arg("pCell"))
        .def("Update",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::Update,
            " ", py::arg("hasHadBirthsOrDeaths") = true)
        .def("TessellateIfNeeded",
            (void(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::TessellateIfNeeded,
            " ")
        .def("DivideLongSprings",
            (void(MeshBasedCellPopulation_2_2::*)(double)) &MeshBasedCellPopulation_2_2::DivideLongSprings,
            " ", py::arg("springDivisionThreshold"))
        .def("GetNode",
            (::Node<2> *(MeshBasedCellPopulation_2_2::*)(unsigned int)) &MeshBasedCellPopulation_2_2::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNumNodes",
            (unsigned int(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetNumNodes,
            " ")
        .def("WriteVtkResultsToFile",
            (void(MeshBasedCellPopulation_2_2::*)(::std::string const &)) &MeshBasedCellPopulation_2_2::WriteVtkResultsToFile,
            " ", py::arg("rDirectory"))
        .def("GetVolumeOfCell",
            (double(MeshBasedCellPopulation_2_2::*)(::CellPtr)) &MeshBasedCellPopulation_2_2::GetVolumeOfCell,
            " ", py::arg("pCell"))
        .def("CreateVoronoiTessellation",
            (void(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::CreateVoronoiTessellation,
            " ")
        .def("GetVolumeOfVoronoiElement",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int)) &MeshBasedCellPopulation_2_2::GetVolumeOfVoronoiElement,
            " ", py::arg("index"))
        .def("GetSurfaceAreaOfVoronoiElement",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int)) &MeshBasedCellPopulation_2_2::GetSurfaceAreaOfVoronoiElement,
            " ", py::arg("index"))
        .def("GetVoronoiEdgeLength",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int, unsigned int)) &MeshBasedCellPopulation_2_2::GetVoronoiEdgeLength,
            " ", py::arg("index1"), py::arg("index2"))
        .def("GetWidth",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int const &)) &MeshBasedCellPopulation_2_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("WriteDataToVisualizerSetupFile",
            (void(MeshBasedCellPopulation_2_2::*)(::out_stream &)) &MeshBasedCellPopulation_2_2::WriteDataToVisualizerSetupFile,
            " ", py::arg("pVizSetupFile"))
        .def("SpringsBegin",
            (::MeshBasedCellPopulation<2>::SpringIterator(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::SpringsBegin,
            " ")
        .def("SpringsEnd",
            (::MeshBasedCellPopulation<2>::SpringIterator(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::SpringsEnd,
            " ")
        .def("CheckCellPointers",
            (void(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::CheckCellPointers,
            " ")
        .def("GetAreaBasedDampingConstantParameter",
            (double(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetAreaBasedDampingConstantParameter,
            " ")
        .def("SetAreaBasedDampingConstantParameter",
            (void(MeshBasedCellPopulation_2_2::*)(double)) &MeshBasedCellPopulation_2_2::SetAreaBasedDampingConstantParameter,
            " ", py::arg("areaBasedDampingConstantParameter"))
        .def("OutputCellPopulationParameters",
            (void(MeshBasedCellPopulation_2_2::*)(::out_stream &)) &MeshBasedCellPopulation_2_2::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("SetWriteVtkAsPoints",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::SetWriteVtkAsPoints,
            " ", py::arg("writeVtkAsPoints"))
        .def("GetWriteVtkAsPoints",
            (bool(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetWriteVtkAsPoints,
            " ")
        .def("SetBoundVoronoiTessellation",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::SetBoundVoronoiTessellation,
            " ", py::arg("boundVoronoiTessellation"))
        .def("GetBoundVoronoiTessellation",
            (bool(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetBoundVoronoiTessellation,
            " ")
        .def("SetScaleBoundByEdgeLength",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::SetScaleBoundByEdgeLength,
            " ", py::arg("scaleBoundByEdgeLength"))
        .def("GetScaleBoundByEdgeLength",
            (bool(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetScaleBoundByEdgeLength,
            " ")
        .def("SetBoundedVoroniTesselationLengthCutoff",
            (void(MeshBasedCellPopulation_2_2::*)(double)) &MeshBasedCellPopulation_2_2::SetBoundedVoroniTesselationLengthCutoff,
            " ", py::arg("boundedVoroniTesselationLengthCutoff"))
        .def("GetBoundedVoroniTesselationLengthCutoff",
            (double(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetBoundedVoroniTesselationLengthCutoff,
            " ")
        .def("SetOffsetNewBoundaryNodes",
            (void(MeshBasedCellPopulation_2_2::*)(bool)) &MeshBasedCellPopulation_2_2::SetOffsetNewBoundaryNodes,
            " ", py::arg("offsetNewBoundaryNodes"))
        .def("GetOffsetNewBoundaryNodes",
            (bool(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::GetOffsetNewBoundaryNodes,
            " ")
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(MeshBasedCellPopulation_2_2::*)(unsigned int)) &MeshBasedCellPopulation_2_2::GetNeighbouringNodeIndices,
            " ", py::arg("index"))
        .def("CalculateRestLengths",
            (void(MeshBasedCellPopulation_2_2::*)()) &MeshBasedCellPopulation_2_2::CalculateRestLengths,
            " ")
        .def("GetRestLength",
            (double(MeshBasedCellPopulation_2_2::*)(unsigned int, unsigned int)) &MeshBasedCellPopulation_2_2::GetRestLength,
            " ", py::arg("indexA"), py::arg("indexB"))
        .def("SetRestLength",
            (void(MeshBasedCellPopulation_2_2::*)(unsigned int, unsigned int, double)) &MeshBasedCellPopulation_2_2::SetRestLength,
            " ", py::arg("indexA"), py::arg("indexB"), py::arg("restLength"))
        .def("AddPopulationWriterVoronoiDataWriter", &MeshBasedCellPopulation_2_2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &MeshBasedCellPopulation_2_2::AddCellWriter<CellLabelWriter>)
    ;
}
