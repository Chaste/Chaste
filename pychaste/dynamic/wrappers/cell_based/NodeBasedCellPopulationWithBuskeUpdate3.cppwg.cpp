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
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"

#include "NodeBasedCellPopulationWithBuskeUpdate3.cppwg.hpp"

namespace py = pybind11;
typedef NodeBasedCellPopulationWithBuskeUpdate<3 > NodeBasedCellPopulationWithBuskeUpdate3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NodeBasedCellPopulationWithBuskeUpdate3_Overrides : public NodeBasedCellPopulationWithBuskeUpdate3{
    public:
    using NodeBasedCellPopulationWithBuskeUpdate3::NodeBasedCellPopulationWithBuskeUpdate;
    void UpdateNodeLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithBuskeUpdate3,
            UpdateNodeLocations,
                    dt);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NodeBasedCellPopulationWithBuskeUpdate3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }

};
void register_NodeBasedCellPopulationWithBuskeUpdate3_class(py::module &m){
py::class_<NodeBasedCellPopulationWithBuskeUpdate3 , NodeBasedCellPopulationWithBuskeUpdate3_Overrides , boost::shared_ptr<NodeBasedCellPopulationWithBuskeUpdate3 >  , NodeBasedCellPopulation<3>  >(m, "NodeBasedCellPopulationWithBuskeUpdate3")
        .def(py::init<::NodesOnlyMesh<3> &, ::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, bool >(), py::arg("rMesh"), py::arg("rCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("deleteMesh") = false)
        .def(py::init<::NodesOnlyMesh<3> & >(), py::arg("rMesh"))
        .def(
            "UpdateNodeLocations",
            (void(NodeBasedCellPopulationWithBuskeUpdate3::*)(double)) &NodeBasedCellPopulationWithBuskeUpdate3::UpdateNodeLocations,
            " " , py::arg("dt") )
        .def(
            "OutputCellPopulationParameters",
            (void(NodeBasedCellPopulationWithBuskeUpdate3::*)(::out_stream &)) &NodeBasedCellPopulationWithBuskeUpdate3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def("AddPopulationWriterVoronoiDataWriter", &NodeBasedCellPopulationWithBuskeUpdate3::AddPopulationWriter<VoronoiDataWriter>)
        .def("AddCellWriterCellLabelWriter", &NodeBasedCellPopulationWithBuskeUpdate3::AddCellWriter<CellLabelWriter>)
    ;
}
