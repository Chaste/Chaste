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
#include "VertexBasedCellPopulation.hpp"

#include "VertexBasedCellPopulation2.cppwg.hpp"

namespace py = pybind11;
typedef VertexBasedCellPopulation<2 > VertexBasedCellPopulation2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::Node<2> * _Node_lt_2_gt_Ptr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::set<std::pair<unsigned int, unsigned int>> _std_set_lt_std_pair_lt_unsignedint_unsignedint_gt__gt_;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::TetrahedralMesh<2, 2> * _TetrahedralMesh_lt_2_2_gt_Ptr;

class VertexBasedCellPopulation2_Overrides : public VertexBasedCellPopulation2{
    public:
    using VertexBasedCellPopulation2::VertexBasedCellPopulation;
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            VertexBasedCellPopulation2,
            GetDampingConstant,
                    nodeIndex);
    }
    unsigned int GetNumNodes() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexBasedCellPopulation2,
            GetNumNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            VertexBasedCellPopulation2,
            GetLocationOfCellCentre,
                    pCell);
    }
    ::Node<2> * GetNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _Node_lt_2_gt_Ptr,
            VertexBasedCellPopulation2,
            GetNode,
                    index);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            VertexBasedCellPopulation2,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    ::std::set<std::pair<unsigned int, unsigned int>> GetNeighbouringEdgeIndices(::CellPtr pCell, unsigned int edgeLocalIndex) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_std_pair_lt_unsignedint_unsignedint_gt__gt_,
            VertexBasedCellPopulation2,
            GetNeighbouringEdgeIndices,
                    pCell,
        edgeLocalIndex);
    }
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexBasedCellPopulation2,
            AddNode,
                    pNewNode);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 2> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> & rNewLocation) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            SetNode,
                    index,
        rNewLocation);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            VertexBasedCellPopulation2,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    unsigned int RemoveDeadCells() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            VertexBasedCellPopulation2,
            RemoveDeadCells,
            );
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            VertexBasedCellPopulation2,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    void Update(bool hasHadBirthsOrDeaths) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            Update,
                    hasHadBirthsOrDeaths);
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            OpenWritersFiles,
                    rOutputFileHandler);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>> pPopulationWriter) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            AcceptPopulationWriter,
                    pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>> pPopulationCountWriter) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            AcceptPopulationCountWriter,
                    pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>> pPopulationEventWriter) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            AcceptPopulationEventWriter,
                    pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<2, 2>> pCellWriter, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            AcceptCellWriter,
                    pCellWriter,
        pCell);
    }
    double GetVolumeOfCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            VertexBasedCellPopulation2,
            GetVolumeOfCell,
                    pCell);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetWidth(unsigned int const & rDimension) override {
        PYBIND11_OVERRIDE(
            double,
            VertexBasedCellPopulation2,
            GetWidth,
                    rDimension);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            VertexBasedCellPopulation2,
            GetNeighbouringNodeIndices,
                    index);
    }
    bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned int pdeNodeIndex) override {
        PYBIND11_OVERRIDE(
            bool,
            VertexBasedCellPopulation2,
            IsPdeNodeAssociatedWithNonApoptoticCell,
                    pdeNodeIndex);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            VertexBasedCellPopulation2,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            VertexBasedCellPopulation2,
            GetDefaultTimeStep,
            );
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            WriteDataToVisualizerSetupFile,
                    pVizSetupFile);
    }
    void SimulationSetupHook(::AbstractCellBasedSimulation<2, 2> * pSimulation) override {
        PYBIND11_OVERRIDE(
            void,
            VertexBasedCellPopulation2,
            SimulationSetupHook,
                    pSimulation);
    }

};
void register_VertexBasedCellPopulation2_class(py::module &m){
py::class_<VertexBasedCellPopulation2 , VertexBasedCellPopulation2_Overrides , boost::shared_ptr<VertexBasedCellPopulation2 >  , AbstractOffLatticeCellPopulation<2>  >(m, "VertexBasedCellPopulation2")
        .def(py::init<::MutableVertexMesh<2, 2> &, ::std::vector<boost::shared_ptr<Cell>> &, bool, bool, ::std::vector<unsigned int> const >(), py::arg("rMesh"), py::arg("rCells"), py::arg("deleteMesh") = false, py::arg("validate") = true, py::arg("locationIndices") = std::vector<unsigned int>())
        .def(py::init<::MutableVertexMesh<2, 2> &, ::VertexBasedPopulationSrn<2> & >(), py::arg("rMesh"), py::arg("rPopSrn"))
        .def(
            "GetDampingConstant",
            (double(VertexBasedCellPopulation2::*)(unsigned int)) &VertexBasedCellPopulation2::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "GetNumElements",
            (unsigned int(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetNumElements,
            " "  )
        .def(
            "GetNumNodes",
            (unsigned int(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetNumNodes,
            " "  )
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 2>(VertexBasedCellPopulation2::*)(::CellPtr)) &VertexBasedCellPopulation2::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "GetNode",
            (::Node<2> *(VertexBasedCellPopulation2::*)(unsigned int)) &VertexBasedCellPopulation2::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(VertexBasedCellPopulation2::*)(::CellPtr)) &VertexBasedCellPopulation2::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "GetNeighbouringEdgeIndices",
            (::std::set<std::pair<unsigned int, unsigned int>>(VertexBasedCellPopulation2::*)(::CellPtr, unsigned int)) &VertexBasedCellPopulation2::GetNeighbouringEdgeIndices,
            " " , py::arg("pCell"), py::arg("edgeLocalIndex") )
        .def(
            "AddNode",
            (unsigned int(VertexBasedCellPopulation2::*)(::Node<2> *)) &VertexBasedCellPopulation2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "CheckForStepSizeException",
            (void(VertexBasedCellPopulation2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2> &, double)) &VertexBasedCellPopulation2::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "SetNode",
            (void(VertexBasedCellPopulation2::*)(unsigned int, ::ChastePoint<2> &)) &VertexBasedCellPopulation2::SetNode,
            " " , py::arg("index"), py::arg("rNewLocation") )
        .def(
            "AddCell",
            (::CellPtr(VertexBasedCellPopulation2::*)(::CellPtr, ::CellPtr)) &VertexBasedCellPopulation2::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "RemoveDeadCells",
            (unsigned int(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::RemoveDeadCells,
            " "  )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(VertexBasedCellPopulation2::*)(::CellPtr)) &VertexBasedCellPopulation2::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "Update",
            (void(VertexBasedCellPopulation2::*)(bool)) &VertexBasedCellPopulation2::Update,
            " " , py::arg("hasHadBirthsOrDeaths") = true )
        .def(
            "OpenWritersFiles",
            (void(VertexBasedCellPopulation2::*)(::OutputFileHandler &)) &VertexBasedCellPopulation2::OpenWritersFiles,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "AcceptPopulationWriter",
            (void(VertexBasedCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationWriter<2, 2>>)) &VertexBasedCellPopulation2::AcceptPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AcceptPopulationCountWriter",
            (void(VertexBasedCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2>>)) &VertexBasedCellPopulation2::AcceptPopulationCountWriter,
            " " , py::arg("pPopulationCountWriter") )
        .def(
            "AcceptPopulationEventWriter",
            (void(VertexBasedCellPopulation2::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2>>)) &VertexBasedCellPopulation2::AcceptPopulationEventWriter,
            " " , py::arg("pPopulationEventWriter") )
        .def(
            "AcceptCellWriter",
            (void(VertexBasedCellPopulation2::*)(::boost::shared_ptr<AbstractCellWriter<2, 2>>, ::CellPtr)) &VertexBasedCellPopulation2::AcceptCellWriter,
            " " , py::arg("pCellWriter"), py::arg("pCell") )
        .def(
            "GetRosetteRankOfCell",
            (unsigned int(VertexBasedCellPopulation2::*)(::CellPtr)) &VertexBasedCellPopulation2::GetRosetteRankOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetVolumeOfCell",
            (double(VertexBasedCellPopulation2::*)(::CellPtr)) &VertexBasedCellPopulation2::GetVolumeOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetLocationsOfT2Swaps",
            (::std::vector<boost::numeric::ublas::c_vector<double, 2>>(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetLocationsOfT2Swaps,
            " "  )
        .def(
            "GetCellIdsOfT2Swaps",
            (::std::vector<unsigned int>(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetCellIdsOfT2Swaps,
            " "  )
        .def(
            "AddLocationOfT2Swap",
            (void(VertexBasedCellPopulation2::*)(::boost::numeric::ublas::c_vector<double, 2>)) &VertexBasedCellPopulation2::AddLocationOfT2Swap,
            " " , py::arg("locationOfT2Swap") )
        .def(
            "AddCellIdOfT2Swap",
            (void(VertexBasedCellPopulation2::*)(unsigned int)) &VertexBasedCellPopulation2::AddCellIdOfT2Swap,
            " " , py::arg("idOfT2Swap") )
        .def(
            "ClearLocationsAndCellIdsOfT2Swaps",
            (void(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::ClearLocationsAndCellIdsOfT2Swaps,
            " "  )
        .def(
            "GetOutputCellRearrangementLocations",
            (bool(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetOutputCellRearrangementLocations,
            " "  )
        .def(
            "SetOutputCellRearrangementLocations",
            (void(VertexBasedCellPopulation2::*)(bool)) &VertexBasedCellPopulation2::SetOutputCellRearrangementLocations,
            " " , py::arg("outputCellRearrangementLocations") )
        .def(
            "OutputCellPopulationParameters",
            (void(VertexBasedCellPopulation2::*)(::out_stream &)) &VertexBasedCellPopulation2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetWidth",
            (double(VertexBasedCellPopulation2::*)(unsigned int const &)) &VertexBasedCellPopulation2::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(VertexBasedCellPopulation2::*)(unsigned int)) &VertexBasedCellPopulation2::GetNeighbouringNodeIndices,
            " " , py::arg("index") )
        .def(
            "IsPdeNodeAssociatedWithNonApoptoticCell",
            (bool(VertexBasedCellPopulation2::*)(unsigned int)) &VertexBasedCellPopulation2::IsPdeNodeAssociatedWithNonApoptoticCell,
            " " , py::arg("pdeNodeIndex") )
        .def(
            "GetCellDataItemAtPdeNode",
            (double(VertexBasedCellPopulation2::*)(unsigned int, ::std::string &, bool, double)) &VertexBasedCellPopulation2::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "GetVertexBasedDivisionRule",
            (::boost::shared_ptr<AbstractVertexBasedDivisionRule<2>>(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetVertexBasedDivisionRule,
            " "  )
        .def(
            "SetVertexBasedDivisionRule",
            (void(VertexBasedCellPopulation2::*)(::boost::shared_ptr<AbstractVertexBasedDivisionRule<2>>)) &VertexBasedCellPopulation2::SetVertexBasedDivisionRule,
            " " , py::arg("pVertexBasedDivisionRule") )
        .def(
            "GetDefaultTimeStep",
            (double(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetDefaultTimeStep,
            " "  )
        .def(
            "WriteDataToVisualizerSetupFile",
            (void(VertexBasedCellPopulation2::*)(::out_stream &)) &VertexBasedCellPopulation2::WriteDataToVisualizerSetupFile,
            " " , py::arg("pVizSetupFile") )
        .def(
            "SimulationSetupHook",
            (void(VertexBasedCellPopulation2::*)(::AbstractCellBasedSimulation<2, 2> *)) &VertexBasedCellPopulation2::SimulationSetupHook,
            " " , py::arg("pSimulation") )
        .def(
            "GetRestrictVertexMovementBoolean",
            (bool(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::GetRestrictVertexMovementBoolean,
            " "  )
        .def(
            "SetRestrictVertexMovementBoolean",
            (void(VertexBasedCellPopulation2::*)(bool)) &VertexBasedCellPopulation2::SetRestrictVertexMovementBoolean,
            " " , py::arg("restrictVertexMovement") )
        .def(
            "rGetVertexBasedPopulationSrn",
            (::VertexBasedPopulationSrn<2> &(VertexBasedCellPopulation2::*)()) &VertexBasedCellPopulation2::rGetVertexBasedPopulationSrn,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetVertexBasedPopulationSrn",
            (::VertexBasedPopulationSrn<2> const &(VertexBasedCellPopulation2::*)() const ) &VertexBasedCellPopulation2::rGetVertexBasedPopulationSrn,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetWriteCellVtkResults",
            (void(VertexBasedCellPopulation2::*)(bool const)) &VertexBasedCellPopulation2::SetWriteCellVtkResults,
            " " , py::arg("new_val") )
        .def(
            "SetWriteEdgeVtkResults",
            (void(VertexBasedCellPopulation2::*)(bool const)) &VertexBasedCellPopulation2::SetWriteEdgeVtkResults,
            " " , py::arg("new_val") )
        .def("AddPopulationWriterVoronoiDataWriter", &VertexBasedCellPopulation2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &VertexBasedCellPopulation2::AddCellWriter<CellLabelWriter>)
    ;
}
