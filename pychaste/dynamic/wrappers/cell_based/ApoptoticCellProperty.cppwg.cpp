#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ApoptoticCellProperty.hpp"

#include "ApoptoticCellProperty.cppwg.hpp"

namespace py = pybind11;
typedef ApoptoticCellProperty ApoptoticCellProperty;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ApoptoticCellProperty_class(py::module &m){
py::class_<ApoptoticCellProperty  , boost::shared_ptr<ApoptoticCellProperty >  , AbstractCellProperty  >(m, "ApoptoticCellProperty")
        .def(py::init<unsigned int >(), py::arg("colour") = 6)
        .def(
            "GetColour",
            (unsigned int(ApoptoticCellProperty::*)() const ) &ApoptoticCellProperty::GetColour,
            " "  )
    ;
}
