#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PdeSimulationTime.hpp"

#include "PdeSimulationTime.cppwg.hpp"

namespace py = pybind11;
typedef PdeSimulationTime PdeSimulationTime;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_PdeSimulationTime_class(py::module &m){
py::class_<PdeSimulationTime  , boost::shared_ptr<PdeSimulationTime >   >(m, "PdeSimulationTime")
        .def(py::init< >())
        .def_static(
            "SetTime",
            (void(*)(double)) &PdeSimulationTime::SetTime,
            " " , py::arg("time") )
        .def_static(
            "GetTime",
            (double(*)()) &PdeSimulationTime::GetTime,
            " "  )
        .def_static(
            "SetPdeTimeStepAndNextTime",
            (void(*)(double, double)) &PdeSimulationTime::SetPdeTimeStepAndNextTime,
            " " , py::arg("timestep"), py::arg("next_time") )
        .def_static(
            "GetPdeTimeStep",
            (double(*)()) &PdeSimulationTime::GetPdeTimeStep,
            " "  )
        .def_static(
            "GetNextTime",
            (double(*)()) &PdeSimulationTime::GetNextTime,
            " "  )
        .def_static(
            "GetPdeTimeStepInverse",
            (double(*)()) &PdeSimulationTime::GetPdeTimeStepInverse,
            " "  )
    ;
}
