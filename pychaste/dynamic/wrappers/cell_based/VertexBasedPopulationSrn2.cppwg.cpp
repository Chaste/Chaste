#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VertexBasedPopulationSrn.hpp"

#include "VertexBasedPopulationSrn2.cppwg.hpp"

namespace py = pybind11;
typedef VertexBasedPopulationSrn<2 > VertexBasedPopulationSrn2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_VertexBasedPopulationSrn2_class(py::module &m){
py::class_<VertexBasedPopulationSrn2  , boost::shared_ptr<VertexBasedPopulationSrn2 >   >(m, "VertexBasedPopulationSrn2")
        .def(py::init< >())
        .def(
            "SetVertexCellPopulation",
            (void(VertexBasedPopulationSrn2::*)(::VertexBasedCellPopulation<2> *)) &VertexBasedPopulationSrn2::SetVertexCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "UpdateSrnAfterBirthOrDeath",
            (void(VertexBasedPopulationSrn2::*)(::VertexElementMap &)) &VertexBasedPopulationSrn2::UpdateSrnAfterBirthOrDeath,
            " " , py::arg("rElementMap") )
        .def(
            "RemapCellSrn",
            (void(VertexBasedPopulationSrn2::*)(::std::vector<boost::shared_ptr<AbstractSrnModel>>, ::CellSrnModel *, ::EdgeRemapInfo const &)) &VertexBasedPopulationSrn2::RemapCellSrn,
            " " , py::arg("parentSrnEdges"), py::arg("pCellSrn"), py::arg("rEdgeChange") )
    ;
}
