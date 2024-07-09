#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchEdgeTrackingModifier.hpp"

#include "DeltaNotchEdgeTrackingModifier2.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeTrackingModifier<2 > DeltaNotchEdgeTrackingModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchEdgeTrackingModifier2_Overrides : public DeltaNotchEdgeTrackingModifier2{
    public:
    using DeltaNotchEdgeTrackingModifier2::DeltaNotchEdgeTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeTrackingModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchEdgeTrackingModifier2_class(py::module &m){
py::class_<DeltaNotchEdgeTrackingModifier2 , DeltaNotchEdgeTrackingModifier2_Overrides , boost::shared_ptr<DeltaNotchEdgeTrackingModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "DeltaNotchEdgeTrackingModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchEdgeTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchEdgeTrackingModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchEdgeTrackingModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &DeltaNotchEdgeTrackingModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchEdgeTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchEdgeTrackingModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchEdgeTrackingModifier2::*)(::out_stream &)) &DeltaNotchEdgeTrackingModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
