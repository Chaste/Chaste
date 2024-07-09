#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeTrackingModifier.hpp"

#include "VolumeTrackingModifier2.cppwg.hpp"

namespace py = pybind11;
typedef VolumeTrackingModifier<2 > VolumeTrackingModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeTrackingModifier2_Overrides : public VolumeTrackingModifier2{
    public:
    using VolumeTrackingModifier2::VolumeTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_VolumeTrackingModifier2_class(py::module &m){
py::class_<VolumeTrackingModifier2 , VolumeTrackingModifier2_Overrides , boost::shared_ptr<VolumeTrackingModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "VolumeTrackingModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(VolumeTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &VolumeTrackingModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(VolumeTrackingModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &VolumeTrackingModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(VolumeTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &VolumeTrackingModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(VolumeTrackingModifier2::*)(::out_stream &)) &VolumeTrackingModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
