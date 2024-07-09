#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VtkSceneModifier.hpp"

#include "VtkSceneModifier2.cppwg.hpp"

namespace py = pybind11;
typedef VtkSceneModifier<2 > VtkSceneModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VtkSceneModifier2_Overrides : public VtkSceneModifier2{
    public:
    using VtkSceneModifier2::VtkSceneModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VtkSceneModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_VtkSceneModifier2_class(py::module &m){
py::class_<VtkSceneModifier2 , VtkSceneModifier2_Overrides , boost::shared_ptr<VtkSceneModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "VtkSceneModifier2")
        .def(py::init< >())
        .def(
            "GetVtkScene",
            (::boost::shared_ptr<VtkScene<2>>(VtkSceneModifier2::*)()) &VtkSceneModifier2::GetVtkScene,
            " "  )
        .def(
            "UpdateAtEndOfTimeStep",
            (void(VtkSceneModifier2::*)(::AbstractCellPopulation<2> &)) &VtkSceneModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(VtkSceneModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &VtkSceneModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetVtkScene",
            (void(VtkSceneModifier2::*)(::boost::shared_ptr<VtkScene<2>>)) &VtkSceneModifier2::SetVtkScene,
            " " , py::arg("pScene") )
        .def(
            "UpdateCellData",
            (void(VtkSceneModifier2::*)(::AbstractCellPopulation<2> &)) &VtkSceneModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(VtkSceneModifier2::*)(::out_stream &)) &VtkSceneModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "SetUpdateFrequency",
            (void(VtkSceneModifier2::*)(unsigned int)) &VtkSceneModifier2::SetUpdateFrequency,
            " " , py::arg("frequency") )
    ;
}
