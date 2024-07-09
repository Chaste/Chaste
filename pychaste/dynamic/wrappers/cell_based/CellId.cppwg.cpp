#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellId.hpp"

#include "CellId.cppwg.hpp"

namespace py = pybind11;
typedef CellId CellId;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellId_class(py::module &m){
py::class_<CellId  , boost::shared_ptr<CellId >  , AbstractCellProperty  >(m, "CellId")
        .def(py::init< >())
        .def(
            "AssignCellId",
            (void(CellId::*)()) &CellId::AssignCellId,
            " "  )
        .def(
            "GetMaxCellId",
            (unsigned int(CellId::*)() const ) &CellId::GetMaxCellId,
            " "  )
        .def(
            "GetCellId",
            (unsigned int(CellId::*)() const ) &CellId::GetCellId,
            " "  )
        .def_static(
            "ResetMaxCellId",
            (void(*)()) &CellId::ResetMaxCellId,
            " "  )
    ;
}
