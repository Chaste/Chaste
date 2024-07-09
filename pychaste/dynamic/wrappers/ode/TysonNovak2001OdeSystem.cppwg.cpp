#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TysonNovak2001OdeSystem.hpp"

#include "TysonNovak2001OdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef TysonNovak2001OdeSystem TysonNovak2001OdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TysonNovak2001OdeSystem_Overrides : public TysonNovak2001OdeSystem{
    public:
    using TysonNovak2001OdeSystem::TysonNovak2001OdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            TysonNovak2001OdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }
    bool CalculateStoppingEvent(double time, ::std::vector<double> const & rY) override {
        PYBIND11_OVERRIDE(
            bool,
            TysonNovak2001OdeSystem,
            CalculateStoppingEvent,
                    time,
        rY);
    }
    double CalculateRootFunction(double time, ::std::vector<double> const & rY) override {
        PYBIND11_OVERRIDE(
            double,
            TysonNovak2001OdeSystem,
            CalculateRootFunction,
                    time,
        rY);
    }

};
void register_TysonNovak2001OdeSystem_class(py::module &m){
py::class_<TysonNovak2001OdeSystem , TysonNovak2001OdeSystem_Overrides , boost::shared_ptr<TysonNovak2001OdeSystem >   >(m, "TysonNovak2001OdeSystem")
        .def(py::init<::std::vector<double> >(), py::arg("stateVariables") = std::vector<double>())
        .def(
            "Init",
            (void(TysonNovak2001OdeSystem::*)()) &TysonNovak2001OdeSystem::Init,
            " "  )
        .def(
            "EvaluateYDerivatives",
            (void(TysonNovak2001OdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &TysonNovak2001OdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
        .def(
            "CalculateStoppingEvent",
            (bool(TysonNovak2001OdeSystem::*)(double, ::std::vector<double> const &)) &TysonNovak2001OdeSystem::CalculateStoppingEvent,
            " " , py::arg("time"), py::arg("rY") )
        .def(
            "CalculateRootFunction",
            (double(TysonNovak2001OdeSystem::*)(double, ::std::vector<double> const &)) &TysonNovak2001OdeSystem::CalculateRootFunction,
            " " , py::arg("time"), py::arg("rY") )
    ;
}
