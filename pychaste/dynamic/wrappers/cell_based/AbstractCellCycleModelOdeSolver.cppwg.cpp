#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellCycleModelOdeSolver.hpp"

#include "AbstractCellCycleModelOdeSolver.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellCycleModelOdeSolver AbstractCellCycleModelOdeSolver;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellCycleModelOdeSolver_Overrides : public AbstractCellCycleModelOdeSolver{
    public:
    using AbstractCellCycleModelOdeSolver::AbstractCellCycleModelOdeSolver;
    bool IsSetUp() override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellCycleModelOdeSolver,
            IsSetUp,
            );
    }
    void Reset() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModelOdeSolver,
            Reset,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModelOdeSolver,
            Initialise,
            );
    }
    bool IsAdaptive() override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellCycleModelOdeSolver,
            IsAdaptive,
            );
    }

};
void register_AbstractCellCycleModelOdeSolver_class(py::module &m){
py::class_<AbstractCellCycleModelOdeSolver , AbstractCellCycleModelOdeSolver_Overrides , boost::shared_ptr<AbstractCellCycleModelOdeSolver >   >(m, "AbstractCellCycleModelOdeSolver")
        .def(py::init< >())
        .def(
            "IsSetUp",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::IsSetUp,
            " "  )
        .def(
            "Reset",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::Reset,
            " "  )
        .def(
            "SolveAndUpdateStateVariable",
            (void(AbstractCellCycleModelOdeSolver::*)(::AbstractOdeSystem *, double, double, double)) &AbstractCellCycleModelOdeSolver::SolveAndUpdateStateVariable,
            " " , py::arg("pAbstractOdeSystem"), py::arg("startTime"), py::arg("endTime"), py::arg("timeStep") )
        .def(
            "Initialise",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::Initialise,
            " "  )
        .def(
            "StoppingEventOccurred",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::StoppingEventOccurred,
            " "  )
        .def(
            "GetStoppingTime",
            (double(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::GetStoppingTime,
            " "  )
        .def(
            "SetSizeOfOdeSystem",
            (void(AbstractCellCycleModelOdeSolver::*)(unsigned int)) &AbstractCellCycleModelOdeSolver::SetSizeOfOdeSystem,
            " " , py::arg("sizeOfOdeSystem") )
        .def(
            "GetSizeOfOdeSystem",
            (unsigned int(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::GetSizeOfOdeSystem,
            " "  )
        .def(
            "CheckForStoppingEvents",
            (void(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::CheckForStoppingEvents,
            " "  )
        .def(
            "SetMaxSteps",
            (void(AbstractCellCycleModelOdeSolver::*)(long int)) &AbstractCellCycleModelOdeSolver::SetMaxSteps,
            " " , py::arg("numSteps") )
        .def(
            "SetTolerances",
            (void(AbstractCellCycleModelOdeSolver::*)(double, double)) &AbstractCellCycleModelOdeSolver::SetTolerances,
            " " , py::arg("relTol") = 1.0E-4, py::arg("absTol") = 9.9999999999999995E-7 )
        .def(
            "IsAdaptive",
            (bool(AbstractCellCycleModelOdeSolver::*)()) &AbstractCellCycleModelOdeSolver::IsAdaptive,
            " "  )
    ;
}
