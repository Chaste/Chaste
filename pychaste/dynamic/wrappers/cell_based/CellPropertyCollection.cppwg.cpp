#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellPropertyCollection.hpp"

#include "CellPropertyCollection.cppwg.hpp"

namespace py = pybind11;
typedef CellPropertyCollection CellPropertyCollection;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellPropertyCollection_class(py::module &m){
py::class_<CellPropertyCollection  , boost::shared_ptr<CellPropertyCollection >   >(m, "CellPropertyCollection")
        .def(py::init< >())
        .def(
            "AddProperty",
            (void(CellPropertyCollection::*)(::boost::shared_ptr<AbstractCellProperty> const &)) &CellPropertyCollection::AddProperty,
            " " , py::arg("rProp") )
        .def(
            "SetCellPropertyRegistry",
            (void(CellPropertyCollection::*)(::CellPropertyRegistry *)) &CellPropertyCollection::SetCellPropertyRegistry,
            " " , py::arg("pRegistry") )
        .def(
            "HasProperty",
            (bool(CellPropertyCollection::*)(::boost::shared_ptr<AbstractCellProperty> const &) const ) &CellPropertyCollection::HasProperty,
            " " , py::arg("rProp") )
        .def(
            "RemoveProperty",
            (void(CellPropertyCollection::*)(::boost::shared_ptr<AbstractCellProperty> const &)) &CellPropertyCollection::RemoveProperty,
            " " , py::arg("rProp") )
        .def(
            "GetSize",
            (unsigned int(CellPropertyCollection::*)() const ) &CellPropertyCollection::GetSize,
            " "  )
        .def(
            "Begin",
            (::CellPropertyCollection::Iterator(CellPropertyCollection::*)()) &CellPropertyCollection::Begin,
            " "  )
        .def(
            "End",
            (::CellPropertyCollection::Iterator(CellPropertyCollection::*)()) &CellPropertyCollection::End,
            " "  )
        .def(
            "GetProperty",
            (::boost::shared_ptr<AbstractCellProperty>(CellPropertyCollection::*)() const ) &CellPropertyCollection::GetProperty,
            " "  )
    ;
}
