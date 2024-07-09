#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchEdgeTrackingModifier.hpp"

#include "DeltaNotchEdgeTrackingModifier3.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeTrackingModifier<3 > DeltaNotchEdgeTrackingModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchEdgeTrackingModifier3_Overrides : public DeltaNotchEdgeTrackingModifier3{
    public:
    using DeltaNotchEdgeTrackingModifier3::DeltaNotchEdgeTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchEdgeTrackingModifier3_class(py::module &m){
py::class_<DeltaNotchEdgeTrackingModifier3 , DeltaNotchEdgeTrackingModifier3_Overrides , boost::shared_ptr<DeltaNotchEdgeTrackingModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "DeltaNotchEdgeTrackingModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchEdgeTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchEdgeTrackingModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchEdgeTrackingModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &DeltaNotchEdgeTrackingModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchEdgeTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchEdgeTrackingModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchEdgeTrackingModifier3::*)(::out_stream &)) &DeltaNotchEdgeTrackingModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
