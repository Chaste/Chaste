#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchOdeSystem.hpp"

#include "DeltaNotchOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchOdeSystem DeltaNotchOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchOdeSystem_Overrides : public DeltaNotchOdeSystem{
    public:
    using DeltaNotchOdeSystem::DeltaNotchOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchOdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }

};
void register_DeltaNotchOdeSystem_class(py::module &m){
py::class_<DeltaNotchOdeSystem , DeltaNotchOdeSystem_Overrides , boost::shared_ptr<DeltaNotchOdeSystem >  , AbstractOdeSystem  >(m, "DeltaNotchOdeSystem")
        .def(py::init<::std::vector<double> >(), py::arg("stateVariables") = std::vector<double>())
        .def(
            "EvaluateYDerivatives",
            (void(DeltaNotchOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &DeltaNotchOdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
    ;
}
