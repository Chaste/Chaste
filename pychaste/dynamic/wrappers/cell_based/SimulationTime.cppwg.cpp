#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SimulationTime.hpp"

#include "SimulationTime.cppwg.hpp"

namespace py = pybind11;
typedef SimulationTime SimulationTime;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_SimulationTime_class(py::module &m){
py::class_<SimulationTime  , boost::shared_ptr<SimulationTime >   >(m, "SimulationTime")
        .def_static(
            "Instance",
            (::SimulationTime *(*)()) &SimulationTime::Instance,
            " "  , py::return_value_policy::reference)
        .def(
            "SetEndTimeAndNumberOfTimeSteps",
            (void(SimulationTime::*)(double, unsigned int)) &SimulationTime::SetEndTimeAndNumberOfTimeSteps,
            " " , py::arg("endTime"), py::arg("totalTimeStepsInSimulation") )
        .def(
            "ResetEndTimeAndNumberOfTimeSteps",
            (void(SimulationTime::*)(double const &, unsigned int const &)) &SimulationTime::ResetEndTimeAndNumberOfTimeSteps,
            " " , py::arg("rEndTime"), py::arg("rNumberOfTimeStepsInThisRun") )
        .def(
            "GetTimeStep",
            (double(SimulationTime::*)() const ) &SimulationTime::GetTimeStep,
            " "  )
        .def(
            "IncrementTimeOneStep",
            (void(SimulationTime::*)()) &SimulationTime::IncrementTimeOneStep,
            " "  )
        .def(
            "GetTimeStepsElapsed",
            (unsigned int(SimulationTime::*)() const ) &SimulationTime::GetTimeStepsElapsed,
            " "  )
        .def(
            "GetTime",
            (double(SimulationTime::*)() const ) &SimulationTime::GetTime,
            " "  )
        .def_static(
            "Destroy",
            (void(*)()) &SimulationTime::Destroy,
            " "  )
        .def(
            "IsStartTimeSetUp",
            (bool(SimulationTime::*)() const ) &SimulationTime::IsStartTimeSetUp,
            " "  )
        .def(
            "IsEndTimeAndNumberOfTimeStepsSetUp",
            (bool(SimulationTime::*)() const ) &SimulationTime::IsEndTimeAndNumberOfTimeStepsSetUp,
            " "  )
        .def(
            "IsFinished",
            (bool(SimulationTime::*)() const ) &SimulationTime::IsFinished,
            " "  )
        .def(
            "SetStartTime",
            (void(SimulationTime::*)(double)) &SimulationTime::SetStartTime,
            " " , py::arg("startTime") )
    ;
}
