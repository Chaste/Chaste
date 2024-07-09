#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeTrackingModifier.hpp"

#include "VolumeTrackingModifier3.cppwg.hpp"

namespace py = pybind11;
typedef VolumeTrackingModifier<3 > VolumeTrackingModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeTrackingModifier3_Overrides : public VolumeTrackingModifier3{
    public:
    using VolumeTrackingModifier3::VolumeTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeTrackingModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_VolumeTrackingModifier3_class(py::module &m){
py::class_<VolumeTrackingModifier3 , VolumeTrackingModifier3_Overrides , boost::shared_ptr<VolumeTrackingModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "VolumeTrackingModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(VolumeTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &VolumeTrackingModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(VolumeTrackingModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &VolumeTrackingModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(VolumeTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &VolumeTrackingModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(VolumeTrackingModifier3::*)(::out_stream &)) &VolumeTrackingModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
