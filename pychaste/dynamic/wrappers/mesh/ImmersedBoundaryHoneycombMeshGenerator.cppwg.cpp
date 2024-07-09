#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"

#include "ImmersedBoundaryHoneycombMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryHoneycombMeshGenerator ImmersedBoundaryHoneycombMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ImmersedBoundaryHoneycombMeshGenerator_class(py::module &m){
py::class_<ImmersedBoundaryHoneycombMeshGenerator  , boost::shared_ptr<ImmersedBoundaryHoneycombMeshGenerator >   >(m, "ImmersedBoundaryHoneycombMeshGenerator")
        .def(py::init<unsigned int, unsigned int, unsigned int, double, double >(), py::arg("numElementsX"), py::arg("numElementsY"), py::arg("numNodesPerEdge"), py::arg("proportionalGap"), py::arg("padding"))
        .def(py::init< >())
        .def(
            "GetMesh",
            (::ImmersedBoundaryMesh<2, 2> *(ImmersedBoundaryHoneycombMeshGenerator::*)()) &ImmersedBoundaryHoneycombMeshGenerator::GetMesh,
            " "  , py::return_value_policy::reference)
        .def(
            "GetUnitHexagon",
            (::std::vector<boost::numeric::ublas::c_vector<double, 2>>(ImmersedBoundaryHoneycombMeshGenerator::*)(unsigned int)) &ImmersedBoundaryHoneycombMeshGenerator::GetUnitHexagon,
            " " , py::arg("numPtsPerSide") )
    ;
}
