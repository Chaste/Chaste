#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchInteriorOdeSystem.hpp"

#include "DeltaNotchInteriorOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchInteriorOdeSystem DeltaNotchInteriorOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchInteriorOdeSystem_Overrides : public DeltaNotchInteriorOdeSystem{
    public:
    using DeltaNotchInteriorOdeSystem::DeltaNotchInteriorOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchInteriorOdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }

};
void register_DeltaNotchInteriorOdeSystem_class(py::module &m){
py::class_<DeltaNotchInteriorOdeSystem , DeltaNotchInteriorOdeSystem_Overrides , boost::shared_ptr<DeltaNotchInteriorOdeSystem >  , AbstractOdeSystem  >(m, "DeltaNotchInteriorOdeSystem")
        .def(py::init<::std::vector<double> >(), py::arg("stateVariables") = std::vector<double>())
        .def(
            "EvaluateYDerivatives",
            (void(DeltaNotchInteriorOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &DeltaNotchInteriorOdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
    ;
}
