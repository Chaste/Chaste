#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellEdgeData.hpp"

#include "CellEdgeData.cppwg.hpp"

namespace py = pybind11;
typedef CellEdgeData CellEdgeData;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellEdgeData_class(py::module &m){
py::class_<CellEdgeData  , boost::shared_ptr<CellEdgeData >  , AbstractCellProperty  >(m, "CellEdgeData")
        .def(py::init< >())
        .def(
            "SetItem",
            (void(CellEdgeData::*)(::std::string const &, ::std::vector<double> const &)) &CellEdgeData::SetItem,
            " " , py::arg("rVariableName"), py::arg("rData") )
        .def(
            "GetItem",
            (::std::vector<double>(CellEdgeData::*)(::std::string const &) const ) &CellEdgeData::GetItem,
            " " , py::arg("rVariableName") )
        .def(
            "GetItemAtIndex",
            (double(CellEdgeData::*)(::std::string const &, unsigned int const)) &CellEdgeData::GetItemAtIndex,
            " " , py::arg("rVariableName"), py::arg("index") )
        .def(
            "GetNumItems",
            (unsigned int(CellEdgeData::*)() const ) &CellEdgeData::GetNumItems,
            " "  )
        .def(
            "GetKeys",
            (::std::vector<std::basic_string<char>>(CellEdgeData::*)() const ) &CellEdgeData::GetKeys,
            " "  )
    ;
}
