#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellProperty.hpp"

#include "AbstractCellProperty.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellProperty AbstractCellProperty;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_AbstractCellProperty_class(py::module &m){
py::class_<AbstractCellProperty  , boost::shared_ptr<AbstractCellProperty >   >(m, "AbstractCellProperty")
        .def(py::init< >())
        .def(
            "IsSame",
            (bool(AbstractCellProperty::*)(::AbstractCellProperty const *) const ) &AbstractCellProperty::IsSame,
            " " , py::arg("pOther") )
        .def(
            "IsSame",
            (bool(AbstractCellProperty::*)(::boost::shared_ptr<const AbstractCellProperty>) const ) &AbstractCellProperty::IsSame,
            " " , py::arg("pOther") )
        .def(
            "IncrementCellCount",
            (void(AbstractCellProperty::*)()) &AbstractCellProperty::IncrementCellCount,
            " "  )
        .def(
            "DecrementCellCount",
            (void(AbstractCellProperty::*)()) &AbstractCellProperty::DecrementCellCount,
            " "  )
        .def(
            "GetCellCount",
            (unsigned int(AbstractCellProperty::*)() const ) &AbstractCellProperty::GetCellCount,
            " "  )
    ;
}
