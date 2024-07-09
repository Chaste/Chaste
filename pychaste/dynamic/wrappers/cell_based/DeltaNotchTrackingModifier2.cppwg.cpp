#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchTrackingModifier.hpp"

#include "DeltaNotchTrackingModifier2.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchTrackingModifier<2 > DeltaNotchTrackingModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchTrackingModifier2_Overrides : public DeltaNotchTrackingModifier2{
    public:
    using DeltaNotchTrackingModifier2::DeltaNotchTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchTrackingModifier2_class(py::module &m){
py::class_<DeltaNotchTrackingModifier2 , DeltaNotchTrackingModifier2_Overrides , boost::shared_ptr<DeltaNotchTrackingModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "DeltaNotchTrackingModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchTrackingModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchTrackingModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &DeltaNotchTrackingModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchTrackingModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchTrackingModifier2::*)(::out_stream &)) &DeltaNotchTrackingModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
