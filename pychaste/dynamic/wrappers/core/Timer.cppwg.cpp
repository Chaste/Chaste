#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Timer.hpp"

#include "Timer.cppwg.hpp"

namespace py = pybind11;
typedef Timer Timer;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_Timer_class(py::module &m){
py::class_<Timer  , boost::shared_ptr<Timer >   >(m, "Timer")
        .def(py::init< >())
        .def_static(
            "Reset",
            (void(*)()) &Timer::Reset,
            " "  )
        .def_static(
            "Print",
            (void(*)(::std::string)) &Timer::Print,
            " " , py::arg("message") )
        .def_static(
            "GetElapsedTime",
            (double(*)()) &Timer::GetElapsedTime,
            " "  )
        .def_static(
            "GetWallTime",
            (double(*)()) &Timer::GetWallTime,
            " "  )
        .def_static(
            "PrintAndReset",
            (void(*)(::std::string)) &Timer::PrintAndReset,
            " " , py::arg("message") )
    ;
}
