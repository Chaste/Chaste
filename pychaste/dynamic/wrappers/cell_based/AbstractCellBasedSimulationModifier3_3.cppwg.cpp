#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

#include "AbstractCellBasedSimulationModifier3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedSimulationModifier<3,3 > AbstractCellBasedSimulationModifier3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellBasedSimulationModifier3_3_Overrides : public AbstractCellBasedSimulationModifier3_3{
    public:
    using AbstractCellBasedSimulationModifier3_3::AbstractCellBasedSimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier3_3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfOutputTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulationModifier3_3,
            UpdateAtEndOfOutputTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier3_3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateAtEndOfSolve(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulationModifier3_3,
            UpdateAtEndOfSolve,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier3_3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractCellBasedSimulationModifier3_3_class(py::module &m){
py::class_<AbstractCellBasedSimulationModifier3_3 , AbstractCellBasedSimulationModifier3_3_Overrides , boost::shared_ptr<AbstractCellBasedSimulationModifier3_3 >   >(m, "AbstractCellBasedSimulationModifier3_3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::AbstractCellPopulation<3> &)) &AbstractCellBasedSimulationModifier3_3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfOutputTimeStep",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::AbstractCellPopulation<3> &)) &AbstractCellBasedSimulationModifier3_3::UpdateAtEndOfOutputTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::AbstractCellPopulation<3> &, ::std::string)) &AbstractCellBasedSimulationModifier3_3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateAtEndOfSolve",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::AbstractCellPopulation<3> &)) &AbstractCellBasedSimulationModifier3_3::UpdateAtEndOfSolve,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierInfo",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::out_stream &)) &AbstractCellBasedSimulationModifier3_3::OutputSimulationModifierInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractCellBasedSimulationModifier3_3::*)(::out_stream &)) &AbstractCellBasedSimulationModifier3_3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
