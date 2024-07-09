#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellAncestor.hpp"

#include "CellAncestor.cppwg.hpp"

namespace py = pybind11;
typedef CellAncestor CellAncestor;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellAncestor_class(py::module &m){
py::class_<CellAncestor  , boost::shared_ptr<CellAncestor >  , AbstractCellProperty  >(m, "CellAncestor")
        .def(py::init<unsigned int >(), py::arg("ancestor") = ::UNSIGNED_UNSET)
        .def(
            "GetAncestor",
            (unsigned int(CellAncestor::*)() const ) &CellAncestor::GetAncestor,
            " "  )
    ;
}
