#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"

#include "Alarcon2004OxygenBasedCellCycleOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef Alarcon2004OxygenBasedCellCycleOdeSystem Alarcon2004OxygenBasedCellCycleOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class Alarcon2004OxygenBasedCellCycleOdeSystem_Overrides : public Alarcon2004OxygenBasedCellCycleOdeSystem{
    public:
    using Alarcon2004OxygenBasedCellCycleOdeSystem::Alarcon2004OxygenBasedCellCycleOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            Alarcon2004OxygenBasedCellCycleOdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }
    bool CalculateStoppingEvent(double time, ::std::vector<double> const & rY) override {
        PYBIND11_OVERRIDE(
            bool,
            Alarcon2004OxygenBasedCellCycleOdeSystem,
            CalculateStoppingEvent,
                    time,
        rY);
    }

};
void register_Alarcon2004OxygenBasedCellCycleOdeSystem_class(py::module &m){
py::class_<Alarcon2004OxygenBasedCellCycleOdeSystem , Alarcon2004OxygenBasedCellCycleOdeSystem_Overrides , boost::shared_ptr<Alarcon2004OxygenBasedCellCycleOdeSystem >  , AbstractOdeSystem  >(m, "Alarcon2004OxygenBasedCellCycleOdeSystem")
        .def(py::init<double, bool, ::std::vector<double> >(), py::arg("oxygenConcentration"), py::arg("isLabelled"), py::arg("stateVariables") = std::vector<double>())
        .def(
            "Init",
            (void(Alarcon2004OxygenBasedCellCycleOdeSystem::*)()) &Alarcon2004OxygenBasedCellCycleOdeSystem::Init,
            " "  )
        .def(
            "EvaluateYDerivatives",
            (void(Alarcon2004OxygenBasedCellCycleOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &Alarcon2004OxygenBasedCellCycleOdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
        .def(
            "CalculateStoppingEvent",
            (bool(Alarcon2004OxygenBasedCellCycleOdeSystem::*)(double, ::std::vector<double> const &)) &Alarcon2004OxygenBasedCellCycleOdeSystem::CalculateStoppingEvent,
            " " , py::arg("time"), py::arg("rY") )
        .def(
            "SetIsLabelled",
            (void(Alarcon2004OxygenBasedCellCycleOdeSystem::*)(bool)) &Alarcon2004OxygenBasedCellCycleOdeSystem::SetIsLabelled,
            " " , py::arg("isLabelled") )
        .def(
            "IsLabelled",
            (bool(Alarcon2004OxygenBasedCellCycleOdeSystem::*)() const ) &Alarcon2004OxygenBasedCellCycleOdeSystem::IsLabelled,
            " "  )
        .def(
            "GetOxygenConcentration",
            (double(Alarcon2004OxygenBasedCellCycleOdeSystem::*)() const ) &Alarcon2004OxygenBasedCellCycleOdeSystem::GetOxygenConcentration,
            " "  )
    ;
}
