#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchEdgeInteriorTrackingModifier.hpp"

#include "DeltaNotchEdgeInteriorTrackingModifier2.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeInteriorTrackingModifier<2 > DeltaNotchEdgeInteriorTrackingModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchEdgeInteriorTrackingModifier2_Overrides : public DeltaNotchEdgeInteriorTrackingModifier2{
    public:
    using DeltaNotchEdgeInteriorTrackingModifier2::DeltaNotchEdgeInteriorTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeInteriorTrackingModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchEdgeInteriorTrackingModifier2_class(py::module &m){
py::class_<DeltaNotchEdgeInteriorTrackingModifier2 , DeltaNotchEdgeInteriorTrackingModifier2_Overrides , boost::shared_ptr<DeltaNotchEdgeInteriorTrackingModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "DeltaNotchEdgeInteriorTrackingModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DeltaNotchEdgeInteriorTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchEdgeInteriorTrackingModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DeltaNotchEdgeInteriorTrackingModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &DeltaNotchEdgeInteriorTrackingModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DeltaNotchEdgeInteriorTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DeltaNotchEdgeInteriorTrackingModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DeltaNotchEdgeInteriorTrackingModifier2::*)(::out_stream &)) &DeltaNotchEdgeInteriorTrackingModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
