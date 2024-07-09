#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "WildTypeCellMutationState.hpp"

#include "WildTypeCellMutationState.cppwg.hpp"

namespace py = pybind11;
typedef WildTypeCellMutationState WildTypeCellMutationState;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_WildTypeCellMutationState_class(py::module &m){
py::class_<WildTypeCellMutationState  , boost::shared_ptr<WildTypeCellMutationState >  , AbstractCellMutationState  >(m, "WildTypeCellMutationState")
        .def(py::init< >())
    ;
}
