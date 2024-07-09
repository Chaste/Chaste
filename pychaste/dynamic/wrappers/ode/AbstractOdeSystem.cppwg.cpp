#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOdeSystem.hpp"

#include "AbstractOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeSystem AbstractOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeSystem_Overrides : public AbstractOdeSystem{
    public:
    using AbstractOdeSystem::AbstractOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }
    bool CalculateStoppingEvent(double time, ::std::vector<double> const & rY) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractOdeSystem,
            CalculateStoppingEvent,
                    time,
        rY);
    }
    double CalculateRootFunction(double time, ::std::vector<double> const & rY) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractOdeSystem,
            CalculateRootFunction,
                    time,
        rY);
    }

};
void register_AbstractOdeSystem_class(py::module &m){
py::class_<AbstractOdeSystem , AbstractOdeSystem_Overrides , boost::shared_ptr<AbstractOdeSystem >   >(m, "AbstractOdeSystem")
        .def(
            "EvaluateYDerivatives",
            (void(AbstractOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &AbstractOdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
        .def(
            "CalculateStoppingEvent",
            (bool(AbstractOdeSystem::*)(double, ::std::vector<double> const &)) &AbstractOdeSystem::CalculateStoppingEvent,
            " " , py::arg("time"), py::arg("rY") )
        .def(
            "CalculateRootFunction",
            (double(AbstractOdeSystem::*)(double, ::std::vector<double> const &)) &AbstractOdeSystem::CalculateRootFunction,
            " " , py::arg("time"), py::arg("rY") )
        .def(
            "GetUseAnalyticJacobian",
            (bool(AbstractOdeSystem::*)()) &AbstractOdeSystem::GetUseAnalyticJacobian,
            " "  )
        .def(
            "rGetConstStateVariables",
            (::std::vector<double> const &(AbstractOdeSystem::*)() const ) &AbstractOdeSystem::rGetConstStateVariables,
            " "  , py::return_value_policy::reference_internal)
    ;
}
