#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearPde.hpp"

#include "AbstractLinearPde3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearPde<3,3 > AbstractLinearPde3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_AbstractLinearPde3_3_class(py::module &m){
py::class_<AbstractLinearPde3_3  , boost::shared_ptr<AbstractLinearPde3_3 >   >(m, "AbstractLinearPde3_3")
        .def(py::init< >())
    ;
}
