#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"

#include "ImmersedBoundaryPalisadeMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryPalisadeMeshGenerator ImmersedBoundaryPalisadeMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ImmersedBoundaryPalisadeMeshGenerator_class(py::module &m){
py::class_<ImmersedBoundaryPalisadeMeshGenerator  , boost::shared_ptr<ImmersedBoundaryPalisadeMeshGenerator >   >(m, "ImmersedBoundaryPalisadeMeshGenerator")
        .def(py::init<unsigned int, unsigned int, double, double, double, bool, bool, bool, unsigned int, double >(), py::arg("numCellsWide"), py::arg("numNodesPerCell") = 100, py::arg("ellipseExponent") = 0.20000000000000001, py::arg("cellAspectRatio") = 2., py::arg("randomYMult") = 0., py::arg("basalLamina") = false, py::arg("apicalLamina") = false, py::arg("leakyLaminas") = false, py::arg("numFluidMeshPoints") = (2147483647 * 2U + 1U), py::arg("absoluteGap") = DOUBLE_UNSET)
        .def(py::init< >())
        .def(
            "GetMesh",
            (::ImmersedBoundaryMesh<2, 2> *(ImmersedBoundaryPalisadeMeshGenerator::*)()) &ImmersedBoundaryPalisadeMeshGenerator::GetMesh,
            " "  , py::return_value_policy::reference)
    ;
}
