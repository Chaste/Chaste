#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonSimulationModifier.hpp"

#include "PythonSimulationModifier3.cppwg.hpp"

namespace py = pybind11;
typedef PythonSimulationModifier<3 > PythonSimulationModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class PythonSimulationModifier3_Overrides : public PythonSimulationModifier3{
    public:
    using PythonSimulationModifier3::PythonSimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateCellData(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier3,
            UpdateCellData,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_PythonSimulationModifier3_class(py::module &m){
py::class_<PythonSimulationModifier3 , PythonSimulationModifier3_Overrides , boost::shared_ptr<PythonSimulationModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "PythonSimulationModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(PythonSimulationModifier3::*)(::AbstractCellPopulation<3> &)) &PythonSimulationModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(PythonSimulationModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &PythonSimulationModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(PythonSimulationModifier3::*)(::AbstractCellPopulation<3> &)) &PythonSimulationModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(PythonSimulationModifier3::*)(::out_stream &)) &PythonSimulationModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
