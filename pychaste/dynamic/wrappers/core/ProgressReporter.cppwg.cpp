#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ProgressReporter.hpp"

#include "ProgressReporter.cppwg.hpp"

namespace py = pybind11;
typedef ProgressReporter ProgressReporter;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_ProgressReporter_class(py::module &m){
py::class_<ProgressReporter  , boost::shared_ptr<ProgressReporter >   >(m, "ProgressReporter")
        .def(py::init<::std::string, double, double >(), py::arg("outputDirectory"), py::arg("startTime"), py::arg("endTime"))
        .def(
            "Update",
            (void(ProgressReporter::*)(double)) &ProgressReporter::Update,
            " " , py::arg("currentTime") )
        .def(
            "PrintFinalising",
            (void(ProgressReporter::*)()) &ProgressReporter::PrintFinalising,
            " "  )
        .def(
            "PrintInitialising",
            (void(ProgressReporter::*)()) &ProgressReporter::PrintInitialising,
            " "  )
    ;
}
