#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellLabel.hpp"

#include "CellLabel.cppwg.hpp"

namespace py = pybind11;
typedef CellLabel CellLabel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_CellLabel_class(py::module &m){
py::class_<CellLabel  , boost::shared_ptr<CellLabel >  , AbstractCellProperty  >(m, "CellLabel")
        .def(py::init<unsigned int >(), py::arg("colour") = 5)
        .def(
            "GetColour",
            (unsigned int(CellLabel::*)() const ) &CellLabel::GetColour,
            " "  )
    ;
}
