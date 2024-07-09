#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VtkSceneModifier.hpp"

#include "VtkSceneModifier3.cppwg.hpp"

namespace py = pybind11;
typedef VtkSceneModifier<3 > VtkSceneModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VtkSceneModifier3_Overrides : public VtkSceneModifier3{
    public:
    using VtkSceneModifier3::VtkSceneModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_VtkSceneModifier3_class(py::module &m){
py::class_<VtkSceneModifier3 , VtkSceneModifier3_Overrides , boost::shared_ptr<VtkSceneModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "VtkSceneModifier3")
        .def(py::init< >())
        .def(
            "GetVtkScene",
            (::boost::shared_ptr<VtkScene<3>>(VtkSceneModifier3::*)()) &VtkSceneModifier3::GetVtkScene,
            " "  )
        .def(
            "UpdateAtEndOfTimeStep",
            (void(VtkSceneModifier3::*)(::AbstractCellPopulation<3> &)) &VtkSceneModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(VtkSceneModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &VtkSceneModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetVtkScene",
            (void(VtkSceneModifier3::*)(::boost::shared_ptr<VtkScene<3>>)) &VtkSceneModifier3::SetVtkScene,
            " " , py::arg("pScene") )
        .def(
            "UpdateCellData",
            (void(VtkSceneModifier3::*)(::AbstractCellPopulation<3> &)) &VtkSceneModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(VtkSceneModifier3::*)(::out_stream &)) &VtkSceneModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetUpdateFrequency",
            (void(VtkSceneModifier3::*)(unsigned int)) &VtkSceneModifier3::SetUpdateFrequency,
            " " , py::arg("frequency") )
    ;
}
