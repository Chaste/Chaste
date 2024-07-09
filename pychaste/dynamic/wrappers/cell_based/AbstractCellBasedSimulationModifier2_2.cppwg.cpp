#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

#include "AbstractCellBasedSimulationModifier2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedSimulationModifier<2,2 > AbstractCellBasedSimulationModifier2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellBasedSimulationModifier2_2_Overrides : public AbstractCellBasedSimulationModifier2_2{
    public:
    using AbstractCellBasedSimulationModifier2_2::AbstractCellBasedSimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier2_2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void UpdateAtEndOfOutputTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulationModifier2_2,
            UpdateAtEndOfOutputTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier2_2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateAtEndOfSolve(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulationModifier2_2,
            UpdateAtEndOfSolve,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulationModifier2_2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractCellBasedSimulationModifier2_2_class(py::module &m){
py::class_<AbstractCellBasedSimulationModifier2_2 , AbstractCellBasedSimulationModifier2_2_Overrides , boost::shared_ptr<AbstractCellBasedSimulationModifier2_2 >   >(m, "AbstractCellBasedSimulationModifier2_2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::AbstractCellPopulation<2> &)) &AbstractCellBasedSimulationModifier2_2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateAtEndOfOutputTimeStep",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::AbstractCellPopulation<2> &)) &AbstractCellBasedSimulationModifier2_2::UpdateAtEndOfOutputTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractCellBasedSimulationModifier2_2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateAtEndOfSolve",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::AbstractCellPopulation<2> &)) &AbstractCellBasedSimulationModifier2_2::UpdateAtEndOfSolve,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierInfo",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::out_stream &)) &AbstractCellBasedSimulationModifier2_2::OutputSimulationModifierInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractCellBasedSimulationModifier2_2::*)(::out_stream &)) &AbstractCellBasedSimulationModifier2_2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
