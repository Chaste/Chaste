#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellPropertyRegistry.hpp"

#include "CellPropertyRegistry.cppwg.hpp"

namespace py = pybind11;
typedef CellPropertyRegistry CellPropertyRegistry;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellPropertyRegistry_class(py::module &m){
py::class_<CellPropertyRegistry  , boost::shared_ptr<CellPropertyRegistry >   >(m, "CellPropertyRegistry")
        .def_static(
            "Instance",
            (::CellPropertyRegistry *(*)()) &CellPropertyRegistry::Instance,
            " "  , py::return_value_policy::reference)
        .def(
            "Clear",
            (void(CellPropertyRegistry::*)()) &CellPropertyRegistry::Clear,
            " "  )
        .def(
            "SpecifyOrdering",
            (void(CellPropertyRegistry::*)(::std::vector<boost::shared_ptr<AbstractCellProperty>> const &)) &CellPropertyRegistry::SpecifyOrdering,
            " " , py::arg("rOrdering") )
        .def(
            "HasOrderingBeenSpecified",
            (bool(CellPropertyRegistry::*)()) &CellPropertyRegistry::HasOrderingBeenSpecified,
            " "  )
    ;
}
