#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonSimulationModifier.hpp"

#include "PythonSimulationModifier2.cppwg.hpp"

namespace py = pybind11;
typedef PythonSimulationModifier<2 > PythonSimulationModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class PythonSimulationModifier2_Overrides : public PythonSimulationModifier2{
    public:
    using PythonSimulationModifier2::PythonSimulationModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateCellData(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier2,
            UpdateCellData,
                    rCellPopulation);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PythonSimulationModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_PythonSimulationModifier2_class(py::module &m){
py::class_<PythonSimulationModifier2 , PythonSimulationModifier2_Overrides , boost::shared_ptr<PythonSimulationModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "PythonSimulationModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(PythonSimulationModifier2::*)(::AbstractCellPopulation<2> &)) &PythonSimulationModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(PythonSimulationModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &PythonSimulationModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(PythonSimulationModifier2::*)(::AbstractCellPopulation<2> &)) &PythonSimulationModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(PythonSimulationModifier2::*)(::out_stream &)) &PythonSimulationModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
