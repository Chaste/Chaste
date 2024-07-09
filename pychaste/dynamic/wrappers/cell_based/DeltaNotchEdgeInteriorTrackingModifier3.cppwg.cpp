#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"

#include "DeltaNotchEdgeInteriorTrackingModifier3.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeInteriorTrackingModifier<3 > DeltaNotchEdgeInteriorTrackingModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchEdgeInteriorTrackingModifier3_Overrides : public DeltaNotchEdgeInteriorTrackingModifier3{
    public:
    using DeltaNotchEdgeInteriorTrackingModifier3::DeltaNotchEdgeInteriorTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchEdgeInteriorTrackingModifier3_class(py::module &m){
py::class_<DeltaNotchEdgeInteriorTrackingModifier3 , DeltaNotchEdgeInteriorTrackingModifier3_Overrides , boost::shared_ptr<DeltaNotchEdgeInteriorTrackingModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "DeltaNotchEdgeInteriorTrackingModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchEdgeInteriorTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchEdgeInteriorTrackingModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchEdgeInteriorTrackingModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &DeltaNotchEdgeInteriorTrackingModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchEdgeInteriorTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DeltaNotchEdgeInteriorTrackingModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchEdgeInteriorTrackingModifier3::*)(::out_stream &)) &DeltaNotchEdgeInteriorTrackingModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
