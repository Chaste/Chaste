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

#include "ImmersedBoundaryCellPopulation2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryCellPopulation<2 > ImmersedBoundaryCellPopulation2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::TetrahedralMesh<2, 2> * _TetrahedralMesh_lt_2_2_gt_Ptr;

class ImmersedBoundaryCellPopulation2_Overrides : public ImmersedBoundaryCellPopulation2{
    public:
    using ImmersedBoundaryCellPopulation2::ImmersedBoundaryCellPopulation;
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation2,
            GetDampingConstant,
                    nodeIndex);
    }
    unsigned int GetNumNodes() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation2,
            GetNumNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            ImmersedBoundaryCellPopulation2,
            GetLocationOfCellCentre,
                    pCell);
    }
    ::Node<2> * GetNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            ImmersedBoundaryCellPopulation2,
            GetNode,
                    index);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            ImmersedBoundaryCellPopulation2,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation2,
            AddNode,
                    pNewNode);
    }
    void UpdateNodeLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            UpdateNodeLocations,
                    dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> & rNewLocation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            SetNode,
                    index,
        rNewLocation);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            ImmersedBoundaryCellPopulation2,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    unsigned int RemoveDeadCells() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ImmersedBoundaryCellPopulation2,
            RemoveDeadCells,
            );
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryCellPopulation2,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            Update,
                    hasHadBirthsOrDeaths);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            OpenWritersFiles,
                    rOutputFileHandler);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>> pPopulationWriter) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            AcceptPopulationWriter,
                    pPopulationWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>> pPopulationEventWriter) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            AcceptPopulationEventWriter,
                    pPopulationEventWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>> pPopulationCountWriter) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            AcceptPopulationCountWriter,
                    pPopulationCountWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<2, 2>> pCellWriter, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            AcceptCellWriter,
                    pCellWriter,
        pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation2,
            GetVolumeOfCell,
                    pCell);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetWidth(unsigned int const & rDimension) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation2,
            GetWidth,
                    rDimension);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            ImmersedBoundaryCellPopulation2,
            GetNeighbouringNodeIndices,
                    index);
    }
    bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned int pdeNodeIndex) override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryCellPopulation2,
            IsPdeNodeAssociatedWithNonApoptoticCell,
                    pdeNodeIndex);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation2,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 2> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryCellPopulation2,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryCellPopulation2,
            GetDefaultTimeStep,
            );
    }

};
void register_ImmersedBoundaryCellPopulation2_class(py::module &m){
py::class_<ImmersedBoundaryCellPopulation2 , ImmersedBoundaryCellPopulation2_Overrides , boost::shared_ptr<ImmersedBoundaryCellPopulation2 >  , AbstractOffLatticeCellPopulation<2>  >(m, "ImmersedBoundaryCellPopulation2")
        .def(py::init<::ImmersedBoundaryMesh<2, 2> &, ::std::vector<boost::shared_ptr<Cell>> &, bool, bool, ::std::vector<unsigned int> const >(), py::arg("rMesh"), py::arg("rCells"), py::arg("deleteMesh") = false, py::arg("validate") = true, py::arg("locationIndices") = std::vector<unsigned int>())
        .def(py::init<::ImmersedBoundaryMesh<2, 2> & >(), py::arg("rMesh"))
        .def(
            "GetDampingConstant",
            (double(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "GetElement",
            (::ImmersedBoundaryElement<2, 2> *(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::GetElement,
            " " , py::arg("elementIndex") , py::return_value_policy::reference)
        .def(
            "GetLamina",
            (::ImmersedBoundaryElement<1, 2> *(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::GetLamina,
            " " , py::arg("laminaIndex") , py::return_value_policy::reference)
        .def(
            "GetNumElements",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::GetNumElements,
            " "  )
        .def(
            "GetNumLaminas",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::GetNumLaminas,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::GetNumNodes,
            " "  )
        .def(
            "SetInteractionDistance",
            (void(ImmersedBoundaryCellPopulation2::*)(double)) &ImmersedBoundaryCellPopulation2::SetInteractionDistance,
            " " , py::arg("newDistance") )
        .def(
            "GetInteractionDistance",
            (double(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::GetInteractionDistance,
            " "  )
        .def(
            "SetReMeshFrequency",
            (void(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::SetReMeshFrequency,
            " " , py::arg("newFrequency") )
        .def(
            "GetReMeshFrequency",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::GetReMeshFrequency,
            " "  )
        .def(
            "SetThrowsStepSizeException",
            (void(ImmersedBoundaryCellPopulation2::*)(bool)) &ImmersedBoundaryCellPopulation2::SetThrowsStepSizeException,
            " " , py::arg("throws") )
        .def(
            "ThrowsStepSizeException",
            (bool(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::ThrowsStepSizeException,
            " "  )
        .def(
            "GetIntrinsicSpacing",
            (double(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::GetIntrinsicSpacing,
            " "  )
        .def(
            "SetCellRearrangementThreshold",
            (void(ImmersedBoundaryCellPopulation2::*)(double)) &ImmersedBoundaryCellPopulation2::SetCellRearrangementThreshold,
            " " , py::arg("newThreshold") )
        .def(
            "GetCellRearrangementThreshold",
            (double(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::GetCellRearrangementThreshold,
            " "  )
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 2>(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "GetNode",
            (::Node<2> *(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "AddNode",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)(::Node<2> *)) &ImmersedBoundaryCellPopulation2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "UpdateNodeLocations",
            (void(ImmersedBoundaryCellPopulation2::*)(double)) &ImmersedBoundaryCellPopulation2::UpdateNodeLocations,
            " " , py::arg("dt") )
        .def(
            "SetNode",
            (void(ImmersedBoundaryCellPopulation2::*)(unsigned int, ::ChastePoint<2> &)) &ImmersedBoundaryCellPopulation2::SetNode,
            " " , py::arg("index"), py::arg("rNewLocation") )
        .def(
            "GetElementCorrespondingToCell",
            (::ImmersedBoundaryElement<2, 2> *(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::GetElementCorrespondingToCell,
            " " , py::arg("pCell") , py::return_value_policy::reference)
        .def(
            "AddCell",
            (::CellPtr(ImmersedBoundaryCellPopulation2::*)(::CellPtr, ::CellPtr)) &ImmersedBoundaryCellPopulation2::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "RemoveDeadCells",
            (unsigned int(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::RemoveDeadCells,
            " "  )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "Update",
            (void(ImmersedBoundaryCellPopulation2::*)(bool)) &ImmersedBoundaryCellPopulation2::Update,
            " " , py::arg("hasHadBirthsOrDeaths") = true )
        .def(
            "OpenWritersFiles",
            (void(ImmersedBoundaryCellPopulation2::*)(::OutputFileHandler &)) &ImmersedBoundaryCellPopulation2::OpenWritersFiles,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "AcceptPopulationWriter",
            (void(ImmersedBoundaryCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>>)) &ImmersedBoundaryCellPopulation2::AcceptPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AcceptPopulationEventWriter",
            (void(ImmersedBoundaryCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>>)) &ImmersedBoundaryCellPopulation2::AcceptPopulationEventWriter,
            " " , py::arg("pPopulationEventWriter") )
        .def(
            "AcceptPopulationCountWriter",
            (void(ImmersedBoundaryCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>>)) &ImmersedBoundaryCellPopulation2::AcceptPopulationCountWriter,
            " " , py::arg("pPopulationCountWriter") )
        .def(
            "AcceptCellWriter",
            (void(ImmersedBoundaryCellPopulation2::*)(::boost::shared_ptr<AbstractCellWriter<2, 2>>, ::CellPtr)) &ImmersedBoundaryCellPopulation2::AcceptCellWriter,
            " " , py::arg("pCellWriter"), py::arg("pCell") )
        .def(
            "GetVolumeOfCell",
            (double(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::GetVolumeOfCell,
            " " , py::arg("pCell") )
        .def(
            "OutputCellPopulationParameters",
            (void(ImmersedBoundaryCellPopulation2::*)(::out_stream &)) &ImmersedBoundaryCellPopulation2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetWidth",
            (double(ImmersedBoundaryCellPopulation2::*)(unsigned int const &)) &ImmersedBoundaryCellPopulation2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::GetNeighbouringNodeIndices,
            " " , py::arg("index") )
        .def(
            "IsPdeNodeAssociatedWithNonApoptoticCell",
            (bool(ImmersedBoundaryCellPopulation2::*)(unsigned int)) &ImmersedBoundaryCellPopulation2::IsPdeNodeAssociatedWithNonApoptoticCell,
            " " , py::arg("pdeNodeIndex") )
        .def(
            "GetCellDataItemAtPdeNode",
            (double(ImmersedBoundaryCellPopulation2::*)(unsigned int, ::std::string &, bool, double)) &ImmersedBoundaryCellPopulation2::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "GetImmersedBoundaryDivisionRule",
            (::boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<2>>(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::GetImmersedBoundaryDivisionRule,
            " "  )
        .def(
            "SetImmersedBoundaryDivisionRule",
            (void(ImmersedBoundaryCellPopulation2::*)(::boost::shared_ptr<AbstractImmersedBoundaryDivisionRule<2>>)) &ImmersedBoundaryCellPopulation2::SetImmersedBoundaryDivisionRule,
            " " , py::arg("pImmersedBoundaryDivisionRule") )
        .def(
            "DoesPopulationHaveActiveSources",
            (bool(ImmersedBoundaryCellPopulation2::*)() const ) &ImmersedBoundaryCellPopulation2::DoesPopulationHaveActiveSources,
            " "  )
        .def(
            "IsCellOnBoundary",
            (bool(ImmersedBoundaryCellPopulation2::*)(::CellPtr)) &ImmersedBoundaryCellPopulation2::IsCellOnBoundary,
            " " , py::arg("pCell") )
        .def(
            "SetIfPopulationHasActiveSources",
            (void(ImmersedBoundaryCellPopulation2::*)(bool)) &ImmersedBoundaryCellPopulation2::SetIfPopulationHasActiveSources,
            " " , py::arg("hasActiveSources") )
        .def(
            "SetOutputNodeRegionToVtk",
            (void(ImmersedBoundaryCellPopulation2::*)(bool)) &ImmersedBoundaryCellPopulation2::SetOutputNodeRegionToVtk,
            " " , py::arg("outputNodeRegionsToVtk") )
        .def(
            "CheckForStepSizeException",
            (void(ImmersedBoundaryCellPopulation2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2> &, double)) &ImmersedBoundaryCellPopulation2::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "GetDefaultTimeStep",
            (double(ImmersedBoundaryCellPopulation2::*)()) &ImmersedBoundaryCellPopulation2::GetDefaultTimeStep,
            " "  )
        .def("AddPopulationWriterVoronoiDataWriter", &ImmersedBoundaryCellPopulation2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &ImmersedBoundaryCellPopulation2::AddCellWriter<CellLabelWriter>)
    ;
}
