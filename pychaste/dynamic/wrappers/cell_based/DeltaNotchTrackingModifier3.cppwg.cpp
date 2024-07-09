#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchTrackingModifier.hpp"

#include "DeltaNotchTrackingModifier3.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchTrackingModifier<3 > DeltaNotchTrackingModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchTrackingModifier3_Overrides : public DeltaNotchTrackingModifier3{
    public:
    using DeltaNotchTrackingModifier3::DeltaNotchTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchTrackingModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchTrackingModifier3_class(py::module &m){
py::class_<DeltaNotchTrackingModifier3 , DeltaNotchTrackingModifier3_Overrides , boost::shared_ptr<DeltaNotchTrackingModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "DeltaNotchTrackingModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchTrackingModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchTrackingModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &DeltaNotchTrackingModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchTrackingModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchTrackingModifier3::*)(::out_stream &)) &DeltaNotchTrackingModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
