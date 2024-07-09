#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "DifferentiatedCellProliferativeType.cppwg.hpp"

namespace py = pybind11;
typedef DifferentiatedCellProliferativeType DifferentiatedCellProliferativeType;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_DifferentiatedCellProliferativeType_class(py::module &m){
py::class_<DifferentiatedCellProliferativeType  , boost::shared_ptr<DifferentiatedCellProliferativeType >  , AbstractCellProliferativeType  >(m, "DifferentiatedCellProliferativeType")
        .def(py::init< >())
    ;
}
