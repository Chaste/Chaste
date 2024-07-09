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

#include "NodeBasedCellPopulationWithParticles2.cppwg.hpp"

namespace py = pybind11;
typedef NodeBasedCellPopulationWithParticles<2 > NodeBasedCellPopulationWithParticles2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::CellPtr _CellPtr;

class NodeBasedCellPopulationWithParticles2_Overrides : public NodeBasedCellPopulationWithParticles2{
    public:
    using NodeBasedCellPopulationWithParticles2::NodeBasedCellPopulationWithParticles;
    void UpdateParticlesAfterReMesh(::NodeMap & rMap) override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles2,
            UpdateParticlesAfterReMesh,
                    rMap);
    }
    bool IsParticle(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            NodeBasedCellPopulationWithParticles2,
            IsParticle,
                    index);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            NodeBasedCellPopulationWithParticles2,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles2,
            WriteVtkResultsToFile,
                    rDirectory);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    void Validate() override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles2,
            Validate,
            );
    }
    void AcceptCellWritersAcrossPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithParticles2,
            AcceptCellWritersAcrossPopulation,
            );
    }

};
void register_NodeBasedCellPopulationWithParticles2_class(py::module &m){
py::class_<NodeBasedCellPopulationWithParticles2 , NodeBasedCellPopulationWithParticles2_Overrides , boost::shared_ptr<NodeBasedCellPopulationWithParticles2 >  , NodeBasedCellPopulation<2>  >(m, "NodeBasedCellPopulationWithParticles2")
        .def(py::init<::NodesOnlyMesh<2> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool >(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("deleteMesh") = false)
        .def(py::init<::NodesOnlyMesh<2> & >(), py::arg("rMesh"))
        .def(
            "UpdateParticlesAfterReMesh",
            (void(NodeBasedCellPopulationWithParticles2::*)(::NodeMap &)) &NodeBasedCellPopulationWithParticles2::UpdateParticlesAfterReMesh,
            " " , py::arg("rMap") )
        .def(
            "IsParticle",
            (bool(NodeBasedCellPopulationWithParticles2::*)(unsigned int)) &NodeBasedCellPopulationWithParticles2::IsParticle,
            " " , py::arg("index") )
        .def(
            "GetParticleIndices",
            (::std::set<unsigned int>(NodeBasedCellPopulationWithParticles2::*)()) &NodeBasedCellPopulationWithParticles2::GetParticleIndices,
            " "  )
        .def(
            "AddCell",
            (::CellPtr(NodeBasedCellPopulationWithParticles2::*)(::CellPtr, ::CellPtr)) &NodeBasedCellPopulationWithParticles2::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") )
        .def(
            "WriteVtkResultsToFile",
            (void(NodeBasedCellPopulationWithParticles2::*)(::std::string const &)) &NodeBasedCellPopulationWithParticles2::WriteVtkResultsToFile,
            " " , py::arg("rDirectory") )
        .def(
            "OutputCellPopulationParameters",
            (void(NodeBasedCellPopulationWithParticles2::*)(::out_stream &)) &NodeBasedCellPopulationWithParticles2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def("AddPopulationWriterVoronoiDataWriter", &NodeBasedCellPopulationWithParticles2::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &NodeBasedCellPopulationWithParticles2::AddCellWriter<CellLabelWriter>)
    ;
}
