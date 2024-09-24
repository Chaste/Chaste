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
#include "PottsBasedCellPopulation.hpp"

#include "PottsBasedCellPopulation_2.cppwg.hpp"

namespace py = pybind11;
typedef PottsBasedCellPopulation<2> PottsBasedCellPopulation_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<2, 2> * _TetrahedralMesh_lt_2_2_gt_Ptr;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;

class PottsBasedCellPopulation_2_Overrides : public PottsBasedCellPopulation_2
{
public:
    using PottsBasedCellPopulation_2::PottsBasedCellPopulation;
    ::Node<2> * GetNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            PottsBasedCellPopulation_2,
            GetNode,
            index);
    }
    unsigned int GetNumNodes() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsBasedCellPopulation_2,
            GetNumNodes,
            );
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            PottsBasedCellPopulation_2,
            GetNeighbouringLocationIndices,
            pCell);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetLocationOfCellCentre(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            PottsBasedCellPopulation_2,
            GetLocationOfCellCentre,
            pCell);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            PottsBasedCellPopulation_2,
            AddCell,
            pNewCell,
            pParentCell);
    }
    unsigned int RemoveDeadCells() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            PottsBasedCellPopulation_2,
            RemoveDeadCells,
            );
    }
    void UpdateCellLocations(double dt) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            UpdateCellLocations,
            dt);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            bool,
            PottsBasedCellPopulation_2,
            IsCellAssociatedWithADeletedLocation,
            pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            Update,
            hasHadBirthsOrDeaths);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            OpenWritersFiles,
            rOutputFileHandler);
    }
    void WriteResultsToFiles(::std::string const & rDirectory) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            WriteResultsToFiles,
            rDirectory);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>> pPopulationWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            AcceptPopulationWriter,
            pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>> pPopulationCountWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            AcceptPopulationCountWriter,
            pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>> pPopulationEventWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            AcceptPopulationEventWriter,
            pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<2, 2>> pCellWriter, ::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            AcceptCellWriter,
            pCellWriter,
            pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation_2,
            GetVolumeOfCell,
            pCell);
    }
    double GetWidth(unsigned int const & rDimension) override
    {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation_2,
            GetWidth,
            rDimension);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<2>> pUpdateRule) override
    {
        PYBIND11_OVERRIDE(
            void,
            PottsBasedCellPopulation_2,
            AddUpdateRule,
            pUpdateRule);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override
    {
        PYBIND11_OVERRIDE(
            double,
            PottsBasedCellPopulation_2,
            GetCellDataItemAtPdeNode,
            pdeNodeIndex,
            rVariableName,
            dirichletBoundaryConditionApplies,
            dirichletBoundaryValue);
    }
};

void register_PottsBasedCellPopulation_2_class(py::module &m)
{
    py::class_<PottsBasedCellPopulation_2, PottsBasedCellPopulation_2_Overrides, boost::shared_ptr<PottsBasedCellPopulation_2>, AbstractOnLatticeCellPopulation<2>>(m, "PottsBasedCellPopulation_2")
        .def(py::init<::PottsMesh<2> &, ::std::vector<boost::shared_ptr<Cell>> &, bool, bool, ::std::vector<unsigned int> const>(), py::arg("rMesh"), py::arg("rCells"), py::arg("deleteMesh") = false, py::arg("validate") = true, py::arg("locationIndices") = std::vector<unsigned int>())
        .def(py::init<::PottsMesh<2> &>(), py::arg("rMesh"))
        .def("GetNumElements",
            (unsigned int(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::GetNumElements,
            " ")
        .def("GetNode",
            (::Node<2> *(PottsBasedCellPopulation_2::*)(unsigned int)) &PottsBasedCellPopulation_2::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNumNodes",
            (unsigned int(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::GetNumNodes,
            " ")
        .def("GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(PottsBasedCellPopulation_2::*)(::CellPtr)) &PottsBasedCellPopulation_2::GetNeighbouringLocationIndices,
            " ", py::arg("pCell"))
        .def("GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 2>(PottsBasedCellPopulation_2::*)(::CellPtr)) &PottsBasedCellPopulation_2::GetLocationOfCellCentre,
            " ", py::arg("pCell"))
        .def("AddCell",
            (::CellPtr(PottsBasedCellPopulation_2::*)(::CellPtr, ::CellPtr)) &PottsBasedCellPopulation_2::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ))
        .def("RemoveDeadCells",
            (unsigned int(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::RemoveDeadCells,
            " ")
        .def("UpdateCellLocations",
            (void(PottsBasedCellPopulation_2::*)(double)) &PottsBasedCellPopulation_2::UpdateCellLocations,
            " ", py::arg("dt"))
        .def("IsCellAssociatedWithADeletedLocation",
            (bool(PottsBasedCellPopulation_2::*)(::CellPtr)) &PottsBasedCellPopulation_2::IsCellAssociatedWithADeletedLocation,
            " ", py::arg("pCell"))
        .def("Update",
            (void(PottsBasedCellPopulation_2::*)(bool)) &PottsBasedCellPopulation_2::Update,
            " ", py::arg("hasHadBirthsOrDeaths") = true)
        .def("OpenWritersFiles",
            (void(PottsBasedCellPopulation_2::*)(::OutputFileHandler &)) &PottsBasedCellPopulation_2::OpenWritersFiles,
            " ", py::arg("rOutputFileHandler"))
        .def("WriteResultsToFiles",
            (void(PottsBasedCellPopulation_2::*)(::std::string const &)) &PottsBasedCellPopulation_2::WriteResultsToFiles,
            " ", py::arg("rDirectory"))
        .def("AcceptPopulationWriter",
            (void(PottsBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>>)) &PottsBasedCellPopulation_2::AcceptPopulationWriter,
            " ", py::arg("pPopulationWriter"))
        .def("AcceptPopulationCountWriter",
            (void(PottsBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>>)) &PottsBasedCellPopulation_2::AcceptPopulationCountWriter,
            " ", py::arg("pPopulationCountWriter"))
        .def("AcceptPopulationEventWriter",
            (void(PottsBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>>)) &PottsBasedCellPopulation_2::AcceptPopulationEventWriter,
            " ", py::arg("pPopulationEventWriter"))
        .def("AcceptCellWriter",
            (void(PottsBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractCellWriter<2, 2>>, ::CellPtr)) &PottsBasedCellPopulation_2::AcceptCellWriter,
            " ", py::arg("pCellWriter"), py::arg("pCell"))
        .def("GetVolumeOfCell",
            (double(PottsBasedCellPopulation_2::*)(::CellPtr)) &PottsBasedCellPopulation_2::GetVolumeOfCell,
            " ", py::arg("pCell"))
        .def("GetWidth",
            (double(PottsBasedCellPopulation_2::*)(unsigned int const &)) &PottsBasedCellPopulation_2::GetWidth,
            " ", py::arg("rDimension"))
        .def("OutputCellPopulationParameters",
            (void(PottsBasedCellPopulation_2::*)(::out_stream &)) &PottsBasedCellPopulation_2::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("SetTemperature",
            (void(PottsBasedCellPopulation_2::*)(double)) &PottsBasedCellPopulation_2::SetTemperature,
            " ", py::arg("temperature"))
        .def("GetTemperature",
            (double(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::GetTemperature,
            " ")
        .def("SetNumSweepsPerTimestep",
            (void(PottsBasedCellPopulation_2::*)(unsigned int)) &PottsBasedCellPopulation_2::SetNumSweepsPerTimestep,
            " ", py::arg("numSweepsPerTimestep"))
        .def("GetNumSweepsPerTimestep",
            (unsigned int(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::GetNumSweepsPerTimestep,
            " ")
        .def("CreateElementTessellation",
            (void(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::CreateElementTessellation,
            " ")
        .def("CreateMutableMesh",
            (void(PottsBasedCellPopulation_2::*)()) &PottsBasedCellPopulation_2::CreateMutableMesh,
            " ")
        .def("AddUpdateRule",
            (void(PottsBasedCellPopulation_2::*)(::boost::shared_ptr<AbstractUpdateRule<2>>)) &PottsBasedCellPopulation_2::AddUpdateRule,
            " ", py::arg("pUpdateRule"))
        .def("GetCellDataItemAtPdeNode",
            (double(PottsBasedCellPopulation_2::*)(unsigned int, ::std::string &, bool, double)) &PottsBasedCellPopulation_2::GetCellDataItemAtPdeNode,
            " ", py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0.)
        .def("AddPopulationWriterVoronoiDataWriter", &PottsBasedCellPopulation_2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &PottsBasedCellPopulation_2::AddCellWriter<CellLabelWriter>)
    ;
}
