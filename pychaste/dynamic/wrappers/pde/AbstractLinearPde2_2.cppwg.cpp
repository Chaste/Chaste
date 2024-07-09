#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractLinearPde.hpp"

#include "AbstractLinearPde2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractLinearPde<2,2 > AbstractLinearPde2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_AbstractLinearPde2_2_class(py::module &m){
py::class_<AbstractLinearPde2_2  , boost::shared_ptr<AbstractLinearPde2_2 >   >(m, "AbstractLinearPde2_2")
        .def(py::init< >())
    ;
}
