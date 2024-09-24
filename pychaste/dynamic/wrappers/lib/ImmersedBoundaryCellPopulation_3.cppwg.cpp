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
#include "ImmersedBoundaryCellPopulation.hpp"

#include "ImmersedBoundaryCellPopulation_3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryCellPopulation<3> ImmersedBoundaryCellPopulation_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;

class ImmersedBoundaryCellPopulation_3_Overrides : public ImmersedBoundaryCellPopulation_3
{
public:
    using ImmersedBoundaryCellPopulation_3::ImmersedBoundaryCellPopulation;
    double GetDampingConstant(unsigned int nodeIndex) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation_3,
            GetDampingConstant,
            nodeIndex);
    }
    unsigned int GetNumNodes() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation_3,
            GetNumNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetLocationOfCellCentre(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            ImmersedBoundaryCellPopulation_3,
            GetLocationOfCellCentre,
            pCell);
    }
    ::Node<3> * GetNode(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _Node_lt_3_gt_Ptr,
            ImmersedBoundaryCellPopulation_3,
            GetNode,
            index);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            ImmersedBoundaryCellPopulation_3,
            GetNeighbouringLocationIndices,
            pCell);
    }
    unsigned int AddNode(::Node<3> * pNewNode) override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation_3,
            AddNode,
            pNewNode);
    }
    void UpdateNodeLocations(double dt) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            UpdateNodeLocations,
            dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<3> & rNewLocation) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            SetNode,
            index,
            rNewLocation);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override
    {
        PYBIND11_OVERRIDE(
            _CellPtr,
            ImmersedBoundaryCellPopulation_3,
            AddCell,
            pNewCell,
            pParentCell);
    }
    unsigned int RemoveDeadCells() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation_3,
            RemoveDeadCells,
            );
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryCellPopulation_3,
            IsCellAssociatedWithADeletedLocation,
            pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            Update,
            hasHadBirthsOrDeaths);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            OpenWritersFiles,
            rOutputFileHandler);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>> pPopulationWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            AcceptPopulationWriter,
            pPopulationWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>> pPopulationEventWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            AcceptPopulationEventWriter,
            pPopulationEventWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>> pPopulationCountWriter) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            AcceptPopulationCountWriter,
            pPopulationCountWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<3, 3>> pCellWriter, ::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            AcceptCellWriter,
            pCellWriter,
            pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation_3,
            GetVolumeOfCell,
            pCell);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            OutputCellPopulationParameters,
            rParamsFile);
    }
    double GetWidth(unsigned int const & rDimension) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation_3,
            GetWidth,
            rDimension);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override
    {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            ImmersedBoundaryCellPopulation_3,
            GetNeighbouringNodeIndices,
            index);
    }
    bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned int pdeNodeIndex) override
    {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryCellPopulation_3,
            IsPdeNodeAssociatedWithNonApoptoticCell,
            pdeNodeIndex);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation_3,
            GetCellDataItemAtPdeNode,
            pdeNodeIndex,
            rVariableName,
            dirichletBoundaryConditionApplies,
            dirichletBoundaryValue);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 3> & rDisplacement, double dt) override
    {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation_3,
            CheckForStepSizeException,
            nodeIndex,
            rDisplacement,
            dt);
    }
    double GetDefaultTimeStep() override
    {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation_3,
            GetDefaultTimeStep,
            );
    }
};

void register_ImmersedBoundaryCellPopulation_3_class(py::module &m)
{
    py::class_<ImmersedBoundaryCellPopulation_3, ImmersedBoundaryCellPopulation_3_Overrides, boost::shared_ptr<ImmersedBoundaryCellPopulation_3>, AbstractOffLatticeCellPopulation<3>>(m, "ImmersedBoundaryCellPopulation_3")
        .def(py::init<::ImmersedBoundaryMesh<3, 3> &, ::std::vector<boost::shared_ptr<Cell>> &, bool, bool, ::std::vector<unsigned int> const>(), py::arg("rMesh"), py::arg("rCells"), py::arg("deleteMesh") = false, py::arg("validate") = true, py::arg("locationIndices") = std::vector<unsigned int>())
        .def(py::init<::ImmersedBoundaryMesh<3, 3> &>(), py::arg("rMesh"))
        .def("GetDampingConstant",
            (double(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::GetDampingConstant,
            " ", py::arg("nodeIndex"))
        .def("GetElement",
            (::ImmersedBoundaryElement<3, 3> *(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::GetElement,
            " ", py::arg("elementIndex"), py::return_value_policy::reference)
        .def("GetLamina",
            (::ImmersedBoundaryElement<2, 3> *(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::GetLamina,
            " ", py::arg("laminaIndex"), py::return_value_policy::reference)
        .def("GetNumElements",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::GetNumElements,
            " ")
        .def("GetNumLaminas",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::GetNumLaminas,
            " ")
        .def("GetNumNodes",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::GetNumNodes,
            " ")
        .def("SetInteractionDistance",
            (void(ImmersedBoundaryCellPopulation_3::*)(double)) &ImmersedBoundaryCellPopulation_3::SetInteractionDistance,
            " ", py::arg("newDistance"))
        .def("GetInteractionDistance",
            (double(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::GetInteractionDistance,
            " ")
        .def("SetReMeshFrequency",
            (void(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::SetReMeshFrequency,
            " ", py::arg("newFrequency"))
        .def("GetReMeshFrequency",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::GetReMeshFrequency,
            " ")
        .def("SetThrowsStepSizeException",
            (void(ImmersedBoundaryCellPopulation_3::*)(bool)) &ImmersedBoundaryCellPopulation_3::SetThrowsStepSizeException,
            " ", py::arg("throws"))
        .def("ThrowsStepSizeException",
            (bool(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::ThrowsStepSizeException,
            " ")
        .def("GetIntrinsicSpacing",
            (double(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::GetIntrinsicSpacing,
            " ")
        .def("SetCellRearrangementThreshold",
            (void(ImmersedBoundaryCellPopulation_3::*)(double)) &ImmersedBoundaryCellPopulation_3::SetCellRearrangementThreshold,
            " ", py::arg("newThreshold"))
        .def("GetCellRearrangementThreshold",
            (double(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::GetCellRearrangementThreshold,
            " ")
        .def("GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 3>(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::GetLocationOfCellCentre,
            " ", py::arg("pCell"))
        .def("GetNode",
            (::Node<3> *(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::GetNode,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::GetNeighbouringLocationIndices,
            " ", py::arg("pCell"))
        .def("AddNode",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)(::Node<3> *)) &ImmersedBoundaryCellPopulation_3::AddNode,
            " ", py::arg("pNewNode"))
        .def("UpdateNodeLocations",
            (void(ImmersedBoundaryCellPopulation_3::*)(double)) &ImmersedBoundaryCellPopulation_3::UpdateNodeLocations,
            " ", py::arg("dt"))
        .def("SetNode",
            (void(ImmersedBoundaryCellPopulation_3::*)(unsigned int, ::ChastePoint<3> &)) &ImmersedBoundaryCellPopulation_3::SetNode,
            " ", py::arg("index"), py::arg("rNewLocation"))
        .def("GetElementCorrespondingToCell",
            (::ImmersedBoundaryElement<3, 3> *(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::GetElementCorrespondingToCell,
            " ", py::arg("pCell"), py::return_value_policy::reference)
        .def("AddCell",
            (::CellPtr(ImmersedBoundaryCellPopulation_3::*)(::CellPtr, ::CellPtr)) &ImmersedBoundaryCellPopulation_3::AddCell,
            " ", py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ))
        .def("RemoveDeadCells",
            (unsigned int(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::RemoveDeadCells,
            " ")
        .def("IsCellAssociatedWithADeletedLocation",
            (bool(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::IsCellAssociatedWithADeletedLocation,
            " ", py::arg("pCell"))
        .def("Update",
            (void(ImmersedBoundaryCellPopulation_3::*)(bool)) &ImmersedBoundaryCellPopulation_3::Update,
            " ", py::arg("hasHadBirthsOrDeaths") = true)
        .def("OpenWritersFiles",
            (void(ImmersedBoundaryCellPopulation_3::*)(::OutputFileHandler &)) &ImmersedBoundaryCellPopulation_3::OpenWritersFiles,
            " ", py::arg("rOutputFileHandler"))
        .def("AcceptPopulationWriter",
            (void(ImmersedBoundaryCellPopulation_3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &ImmersedBoundaryCellPopulation_3::AcceptPopulationWriter,
            " ", py::arg("pPopulationWriter"))
        .def("AcceptPopulationEventWriter",
            (void(ImmersedBoundaryCellPopulation_3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &ImmersedBoundaryCellPopulation_3::AcceptPopulationEventWriter,
            " ", py::arg("pPopulationEventWriter"))
        .def("AcceptPopulationCountWriter",
            (void(ImmersedBoundaryCellPopulation_3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &ImmersedBoundaryCellPopulation_3::AcceptPopulationCountWriter,
            " ", py::arg("pPopulationCountWriter"))
        .def("AcceptCellWriter",
            (void(ImmersedBoundaryCellPopulation_3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>, ::CellPtr)) &ImmersedBoundaryCellPopulation_3::AcceptCellWriter,
            " ", py::arg("pCellWriter"), py::arg("pCell"))
        .def("GetVolumeOfCell",
            (double(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::GetVolumeOfCell,
            " ", py::arg("pCell"))
        .def("OutputCellPopulationParameters",
            (void(ImmersedBoundaryCellPopulation_3::*)(::out_stream &)) &ImmersedBoundaryCellPopulation_3::OutputCellPopulationParameters,
            " ", py::arg("rParamsFile"))
        .def("GetWidth",
            (double(ImmersedBoundaryCellPopulation_3::*)(unsigned int const &)) &ImmersedBoundaryCellPopulation_3::GetWidth,
            " ", py::arg("rDimension"))
        .def("GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::GetNeighbouringNodeIndices,
            " ", py::arg("index"))
        .def("IsPdeNodeAssociatedWithNonApoptoticCell",
            (bool(ImmersedBoundaryCellPopulation_3::*)(unsigned int)) &ImmersedBoundaryCellPopulation_3::IsPdeNodeAssociatedWithNonApoptoticCell,
            " ", py::arg("pdeNodeIndex"))
        .def("GetCellDataItemAtPdeNode",
            (double(ImmersedBoundaryCellPopulation_3::*)(unsigned int, ::std::string &, bool, double)) &ImmersedBoundaryCellPopulation_3::GetCellDataItemAtPdeNode,
            " ", py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0.)
        .def("GetImmersedBoundaryDivisionRule",
            (::boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<3>>(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::GetImmersedBoundaryDivisionRule,
            " ")
        .def("SetImmersedBoundaryDivisionRule",
            (void(ImmersedBoundaryCellPopulation_3::*)(::boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<3>>)) &ImmersedBoundaryCellPopulation_3::SetImmersedBoundaryDivisionRule,
            " ", py::arg("pImmersedBoundaryDivisionRule"))
        .def("DoesPopulationHaveActiveSources",
            (bool(ImmersedBoundaryCellPopulation_3::*)() const) &ImmersedBoundaryCellPopulation_3::DoesPopulationHaveActiveSources,
            " ")
        .def("IsCellOnBoundary",
            (bool(ImmersedBoundaryCellPopulation_3::*)(::CellPtr)) &ImmersedBoundaryCellPopulation_3::IsCellOnBoundary,
            " ", py::arg("pCell"))
        .def("SetIfPopulationHasActiveSources",
            (void(ImmersedBoundaryCellPopulation_3::*)(bool)) &ImmersedBoundaryCellPopulation_3::SetIfPopulationHasActiveSources,
            " ", py::arg("hasActiveSources"))
        .def("SetOutputNodeRegionToVtk",
            (void(ImmersedBoundaryCellPopulation_3::*)(bool)) &ImmersedBoundaryCellPopulation_3::SetOutputNodeRegionToVtk,
            " ", py::arg("outputNodeRegionsToVtk"))
        .def("CheckForStepSizeException",
            (void(ImmersedBoundaryCellPopulation_3::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 3> &, double)) &ImmersedBoundaryCellPopulation_3::CheckForStepSizeException,
            " ", py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt"))
        .def("GetDefaultTimeStep",
            (double(ImmersedBoundaryCellPopulation_3::*)()) &ImmersedBoundaryCellPopulation_3::GetDefaultTimeStep,
            " ")
        .def("AddPopulationWriterVoronoiDataWriter", &ImmersedBoundaryCellPopulation_3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &ImmersedBoundaryCellPopulation_3::AddCellWriter<CellLabelWriter>)
    ;
}
