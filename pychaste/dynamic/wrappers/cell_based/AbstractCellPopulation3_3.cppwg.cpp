#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellBasedSimulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellPopulation.hpp"

#include "AbstractCellPopulation3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellPopulation<3,3 > AbstractCellPopulation3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::TetrahedralMesh<3, 3> * _TetrahedralMesh_lt_3_3_gt_Ptr;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::Node<3> * _Node_lt_3_gt_Ptr;
typedef ::CellPtr _CellPtr;
typedef unsigned int unsignedint;
typedef ::CellPtr _CellPtr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::set<std::pair<unsigned int, unsigned int>> _std_set_lt_std_pair_lt_unsignedint_unsignedint_gt__gt_;

class AbstractCellPopulation3_3_Overrides : public AbstractCellPopulation3_3{
    public:
    using AbstractCellPopulation3_3::AbstractCellPopulation;
    ::TetrahedralMesh<3, 3> * GetTetrahedralMeshForPdeModifier() override {
        PYBIND11_OVERRIDE_PURE(
            _TetrahedralMesh_lt_3_3_gt_Ptr,
            AbstractCellPopulation3_3,
            GetTetrahedralMeshForPdeModifier,
            );
    }
    bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned int pdeNodeIndex) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellPopulation3_3,
            IsPdeNodeAssociatedWithNonApoptoticCell,
                    pdeNodeIndex);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellPopulation3_3,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    unsigned int GetNumNodes() override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractCellPopulation3_3,
            GetNumNodes,
            );
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractCellPopulation3_3,
            GetLocationOfCellCentre,
                    pCell);
    }
    ::Node<3> * GetNode(unsigned int index) override {
        PYBIND11_OVERRIDE_PURE(
            _Node_lt_3_gt_Ptr,
            AbstractCellPopulation3_3,
            GetNode,
                    index);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> & rNewLocation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            SetNode,
                    nodeIndex,
        rNewLocation);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellPopulation3_3,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            WriteDataToVisualizerSetupFile,
                    pVizSetupFile);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE_PURE(
            _CellPtr,
            AbstractCellPopulation3_3,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellPopulation3_3,
            GetDefaultTimeStep,
            );
    }
    unsigned int RemoveDeadCells() override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractCellPopulation3_3,
            RemoveDeadCells,
            );
    }
    void Update(bool hasHadBirthsOrDeaths) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            Update,
                    hasHadBirthsOrDeaths);
    }
    ::CellPtr GetCellUsingLocationIndex(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            AbstractCellPopulation3_3,
            GetCellUsingLocationIndex,
                    index);
    }
    bool IsCellAttachedToLocationIndex(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellPopulation3_3,
            IsCellAttachedToLocationIndex,
                    index);
    }
    void AddCellUsingLocationIndex(unsigned int index, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            AddCellUsingLocationIndex,
                    index,
        pCell);
    }
    void RemoveCellUsingLocationIndex(unsigned int index, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            RemoveCellUsingLocationIndex,
                    index,
        pCell);
    }
    double GetWidth(unsigned int const & rDimension) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellPopulation3_3,
            GetWidth,
                    rDimension);
    }
    double GetVolumeOfCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellPopulation3_3,
            GetVolumeOfCell,
                    pCell);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override {
        PYBIND11_OVERRIDE_PURE(
            _std_set_lt_unsignedint_gt_,
            AbstractCellPopulation3_3,
            GetNeighbouringNodeIndices,
                    index);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE_PURE(
            _std_set_lt_unsignedint_gt_,
            AbstractCellPopulation3_3,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    ::std::set<std::pair<unsigned int, unsigned int>> GetNeighbouringEdgeIndices(::CellPtr pCell, unsigned int pEdgeIndex) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_std_pair_lt_unsignedint_unsignedint_gt__gt_,
            AbstractCellPopulation3_3,
            GetNeighbouringEdgeIndices,
                    pCell,
        pEdgeIndex);
    }
    void UpdateCellProcessLocation() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            UpdateCellProcessLocation,
            );
    }
    void OpenWritersFiles(::OutputFileHandler & rOutputFileHandler) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            OpenWritersFiles,
                    rOutputFileHandler);
    }
    void CloseWritersFiles() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            CloseWritersFiles,
            );
    }
    void WriteResultsToFiles(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            WriteResultsToFiles,
                    rDirectory);
    }
    void AcceptPopulationWriter(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>> pPopulationWriter) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            AcceptPopulationWriter,
                    pPopulationWriter);
    }
    void AcceptPopulationCountWriter(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>> pPopulationCountWriter) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            AcceptPopulationCountWriter,
                    pPopulationCountWriter);
    }
    void AcceptPopulationEventWriter(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>> pPopulationEventWriter) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            AcceptPopulationEventWriter,
                    pPopulationEventWriter);
    }
    void AcceptCellWriter(::boost::shared_ptr<AbstractCellWriter<3, 3>> pCellWriter, ::CellPtr pCell) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            AcceptCellWriter,
                    pCellWriter,
        pCell);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    void SimulationSetupHook(::AbstractCellBasedSimulation<3, 3> * pSimulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            SimulationSetupHook,
                    pSimulation);
    }
    bool IsRoomToDivide(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellPopulation3_3,
            IsRoomToDivide,
                    pCell);
    }
    void Validate() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            Validate,
            );
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulation3_3,
            WriteVtkResultsToFile,
                    rDirectory);
    }
    void AcceptCellWritersAcrossPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulation3_3,
            AcceptCellWritersAcrossPopulation,
            );
    }

};
void register_AbstractCellPopulation3_3_class(py::module &m){
py::class_<AbstractCellPopulation3_3 , AbstractCellPopulation3_3_Overrides , boost::shared_ptr<AbstractCellPopulation3_3 >   >(m, "AbstractCellPopulation3_3")
        .def(py::init<::AbstractMesh<3, 3> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const >(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>())
        .def(
            "InitialiseCells",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::InitialiseCells,
            " "  )
        .def(
            "SetDataOnAllCells",
            (void(AbstractCellPopulation3_3::*)(::std::string const &, double)) &AbstractCellPopulation3_3::SetDataOnAllCells,
            " " , py::arg("rDataName"), py::arg("dataValue") )
        .def(
            "rGetMesh",
            (::AbstractMesh<3, 3> &(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::rGetMesh,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetTetrahedralMeshForPdeModifier",
            (::TetrahedralMesh<3, 3> *(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetTetrahedralMeshForPdeModifier,
            " "  , py::return_value_policy::reference)
        .def(
            "IsPdeNodeAssociatedWithNonApoptoticCell",
            (bool(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::IsPdeNodeAssociatedWithNonApoptoticCell,
            " " , py::arg("pdeNodeIndex") )
        .def(
            "GetCellDataItemAtPdeNode",
            (double(AbstractCellPopulation3_3::*)(unsigned int, ::std::string &, bool, double)) &AbstractCellPopulation3_3::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "GetNumNodes",
            (unsigned int(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetNumNodes,
            " "  )
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "GetNode",
            (::Node<3> *(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::GetNode,
            " " , py::arg("index") , py::return_value_policy::reference)
        .def(
            "SetNode",
            (void(AbstractCellPopulation3_3::*)(unsigned int, ::ChastePoint<3> &)) &AbstractCellPopulation3_3::SetNode,
            " " , py::arg("nodeIndex"), py::arg("rNewLocation") )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "WriteDataToVisualizerSetupFile",
            (void(AbstractCellPopulation3_3::*)(::out_stream &)) &AbstractCellPopulation3_3::WriteDataToVisualizerSetupFile,
            " " , py::arg("pVizSetupFile") )
        .def(
            "AddCell",
            (::CellPtr(AbstractCellPopulation3_3::*)(::CellPtr, ::CellPtr)) &AbstractCellPopulation3_3::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "GetDefaultTimeStep",
            (double(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetDefaultTimeStep,
            " "  )
        .def(
            "RemoveDeadCells",
            (unsigned int(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::RemoveDeadCells,
            " "  )
        .def(
            "Update",
            (void(AbstractCellPopulation3_3::*)(bool)) &AbstractCellPopulation3_3::Update,
            " " , py::arg("hasHadBirthsOrDeaths") = true )
        .def(
            "GetCellMutationStateCount",
            (::std::vector<unsigned int>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCellMutationStateCount,
            " "  )
        .def(
            "GetCellProliferativeTypeCount",
            (::std::vector<unsigned int>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCellProliferativeTypeCount,
            " "  )
        .def(
            "GetCellCyclePhaseCount",
            (::std::vector<unsigned int>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCellCyclePhaseCount,
            " "  )
        .def(
            "GetNumRealCells",
            (unsigned int(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetNumRealCells,
            " "  )
        .def(
            "GetNumAllCells",
            (unsigned int(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetNumAllCells,
            " "  )
        .def(
            "SetCellAncestorsToLocationIndices",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::SetCellAncestorsToLocationIndices,
            " "  )
        .def(
            "GetCellAncestors",
            (::std::set<unsigned int>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCellAncestors,
            " "  )
        .def(
            "GetCellUsingLocationIndex",
            (::CellPtr(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::GetCellUsingLocationIndex,
            " " , py::arg("index") )
        .def(
            "GetCellsUsingLocationIndex",
            (::std::set<boost::shared_ptr<Cell>>(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::GetCellsUsingLocationIndex,
            " " , py::arg("index") )
        .def(
            "IsCellAttachedToLocationIndex",
            (bool(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::IsCellAttachedToLocationIndex,
            " " , py::arg("index") )
        .def(
            "SetCellUsingLocationIndex",
            (void(AbstractCellPopulation3_3::*)(unsigned int, ::CellPtr)) &AbstractCellPopulation3_3::SetCellUsingLocationIndex,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "AddCellUsingLocationIndex",
            (void(AbstractCellPopulation3_3::*)(unsigned int, ::CellPtr)) &AbstractCellPopulation3_3::AddCellUsingLocationIndex,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "RemoveCellUsingLocationIndex",
            (void(AbstractCellPopulation3_3::*)(unsigned int, ::CellPtr)) &AbstractCellPopulation3_3::RemoveCellUsingLocationIndex,
            " " , py::arg("index"), py::arg("pCell") )
        .def(
            "MoveCellInLocationMap",
            (void(AbstractCellPopulation3_3::*)(::CellPtr, unsigned int, unsigned int)) &AbstractCellPopulation3_3::MoveCellInLocationMap,
            " " , py::arg("pCell"), py::arg("old_index"), py::arg("new_index") )
        .def(
            "GetLocationIndexUsingCell",
            (unsigned int(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::GetLocationIndexUsingCell,
            " " , py::arg("pCell") )
        .def(
            "GetCellPropertyRegistry",
            (::boost::shared_ptr<CellPropertyRegistry>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCellPropertyRegistry,
            " "  )
        .def(
            "SetDefaultCellMutationStateAndProliferativeTypeOrdering",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::SetDefaultCellMutationStateAndProliferativeTypeOrdering,
            " "  )
        .def(
            "GetWidth",
            (double(AbstractCellPopulation3_3::*)(unsigned int const &)) &AbstractCellPopulation3_3::GetWidth,
            " " , py::arg("rDimension") )
        .def(
            "GetVolumeOfCell",
            (double(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::GetVolumeOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(AbstractCellPopulation3_3::*)(unsigned int)) &AbstractCellPopulation3_3::GetNeighbouringNodeIndices,
            " " , py::arg("index") )
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "GetNeighbouringEdgeIndices",
            (::std::set<std::pair<unsigned int, unsigned int>>(AbstractCellPopulation3_3::*)(::CellPtr, unsigned int)) &AbstractCellPopulation3_3::GetNeighbouringEdgeIndices,
            " " , py::arg("pCell"), py::arg("pEdgeIndex") )
        .def(
            "GetCentroidOfCellPopulation",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetCentroidOfCellPopulation,
            " "  )
        .def(
            "UpdateCellProcessLocation",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::UpdateCellProcessLocation,
            " "  )
        .def(
            "OpenWritersFiles",
            (void(AbstractCellPopulation3_3::*)(::OutputFileHandler &)) &AbstractCellPopulation3_3::OpenWritersFiles,
            " " , py::arg("rOutputFileHandler") )
        .def(
            "CloseWritersFiles",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::CloseWritersFiles,
            " "  )
        .def(
            "WriteResultsToFiles",
            (void(AbstractCellPopulation3_3::*)(::std::string const &)) &AbstractCellPopulation3_3::WriteResultsToFiles,
            " " , py::arg("rDirectory") )
        .def(
            "AcceptPopulationWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &AbstractCellPopulation3_3::AcceptPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AcceptPopulationCountWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &AbstractCellPopulation3_3::AcceptPopulationCountWriter,
            " " , py::arg("pPopulationCountWriter") )
        .def(
            "AcceptPopulationEventWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &AbstractCellPopulation3_3::AcceptPopulationEventWriter,
            " " , py::arg("pPopulationEventWriter") )
        .def(
            "AcceptCellWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>, ::CellPtr)) &AbstractCellPopulation3_3::AcceptCellWriter,
            " " , py::arg("pCellWriter"), py::arg("pCell") )
        .def(
            "GetDivisionsInformation",
            (::std::vector<std::basic_string<char>>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetDivisionsInformation,
            " "  )
        .def(
            "AddDivisionInformation",
            (void(AbstractCellPopulation3_3::*)(::std::string)) &AbstractCellPopulation3_3::AddDivisionInformation,
            " " , py::arg("divisionInformation") )
        .def(
            "ClearDivisionsInformation",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::ClearDivisionsInformation,
            " "  )
        .def(
            "GetRemovalsInformation",
            (::std::vector<std::basic_string<char>>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetRemovalsInformation,
            " "  )
        .def(
            "AddRemovalInformation",
            (void(AbstractCellPopulation3_3::*)(::std::string)) &AbstractCellPopulation3_3::AddRemovalInformation,
            " " , py::arg("removalInformation") )
        .def(
            "ClearRemovalsInformation",
            (void(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::ClearRemovalsInformation,
            " "  )
        .def(
            "GenerateRemovalInformation",
            (void(AbstractCellPopulation3_3::*)(::CellPtr, ::std::string)) &AbstractCellPopulation3_3::GenerateRemovalInformation,
            " " , py::arg("pCell"), py::arg("killerInfo") )
        .def(
            "KillCell",
            (void(AbstractCellPopulation3_3::*)(::CellPtr, ::std::string)) &AbstractCellPopulation3_3::KillCell,
            " " , py::arg("pCell"), py::arg("killerInfo") )
        .def(
            "StartApoptosisOnCell",
            (void(AbstractCellPopulation3_3::*)(::CellPtr, ::std::string)) &AbstractCellPopulation3_3::StartApoptosisOnCell,
            " " , py::arg("pCell"), py::arg("killerInfo") )
        .def(
            "OutputCellPopulationInfo",
            (void(AbstractCellPopulation3_3::*)(::out_stream &)) &AbstractCellPopulation3_3::OutputCellPopulationInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractCellPopulation3_3::*)(::out_stream &)) &AbstractCellPopulation3_3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SimulationSetupHook",
            (void(AbstractCellPopulation3_3::*)(::AbstractCellBasedSimulation<3, 3> *)) &AbstractCellPopulation3_3::SimulationSetupHook,
            " " , py::arg("pSimulation") )
        .def(
            "GetOutputResultsForChasteVisualizer",
            (bool(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetOutputResultsForChasteVisualizer,
            " "  )
        .def(
            "AddPopulationWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationWriter<3, 3>>)) &AbstractCellPopulation3_3::AddPopulationWriter,
            " " , py::arg("pPopulationWriter") )
        .def(
            "AddCellWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellWriter<3, 3>>)) &AbstractCellPopulation3_3::AddCellWriter,
            " " , py::arg("pCellWriter") )
        .def(
            "AddCellPopulationCountWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationCountWriter<3, 3>>)) &AbstractCellPopulation3_3::AddCellPopulationCountWriter,
            " " , py::arg("pCellPopulationCountWriter") )
        .def(
            "AddCellPopulationEventWriter",
            (void(AbstractCellPopulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationEventWriter<3, 3>>)) &AbstractCellPopulation3_3::AddCellPopulationEventWriter,
            " " , py::arg("pCellPopulationEventWriter") )
        .def(
            "SetOutputResultsForChasteVisualizer",
            (void(AbstractCellPopulation3_3::*)(bool)) &AbstractCellPopulation3_3::SetOutputResultsForChasteVisualizer,
            " " , py::arg("outputResultsForChasteVisualizer") )
        .def(
            "GetSizeOfCellPopulation",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::GetSizeOfCellPopulation,
            " "  )
        .def(
            "IsRoomToDivide",
            (bool(AbstractCellPopulation3_3::*)(::CellPtr)) &AbstractCellPopulation3_3::IsRoomToDivide,
            " " , py::arg("pCell") )
        .def(
            "CreateOrderedPair",
            (::std::pair<unsigned int, unsigned int>(AbstractCellPopulation3_3::*)(unsigned int, unsigned int)) &AbstractCellPopulation3_3::CreateOrderedPair,
            " " , py::arg("index1"), py::arg("index2") )
        .def(
            "Begin",
            (::AbstractCellPopulation<3>::Iterator(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::Begin,
            " "  )
        .def(
            "End",
            (::AbstractCellPopulation<3>::Iterator(AbstractCellPopulation3_3::*)()) &AbstractCellPopulation3_3::End,
            " "  )
    ;
}
