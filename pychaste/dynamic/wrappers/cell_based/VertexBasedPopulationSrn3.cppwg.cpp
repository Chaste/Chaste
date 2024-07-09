#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VertexBasedPopulationSrn.hpp"

#include "VertexBasedPopulationSrn3.cppwg.hpp"

namespace py = pybind11;
typedef VertexBasedPopulationSrn<3 > VertexBasedPopulationSrn3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_VertexBasedPopulationSrn3_class(py::module &m){
py::class_<VertexBasedPopulationSrn3  , boost::shared_ptr<VertexBasedPopulationSrn3 >   >(m, "VertexBasedPopulationSrn3")
        .def(py::init< >())
        .def(
            "SetVertexCellPopulation",
            (void(VertexBasedPopulationSrn3::*)(::VertexBasedCellPopulation<3> *)) &VertexBasedPopulationSrn3::SetVertexCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "UpdateSrnAfterBirthOrDeath",
            (void(VertexBasedPopulationSrn3::*)(::VertexElementMap &)) &VertexBasedPopulationSrn3::UpdateSrnAfterBirthOrDeath,
            " " , py::arg("rElementMap") )
        .def(
            "RemapCellSrn",
            (void(VertexBasedPopulationSrn3::*)(::std::vector<boost::shared_ptr<AbstractSrnModel>>, ::CellSrnModel *, ::EdgeRemapInfo const &)) &VertexBasedPopulationSrn3::RemapCellSrn,
            " " , py::arg("parentSrnEdges"), py::arg("pCellSrn"), py::arg("rEdgeChange") )
    ;
}
