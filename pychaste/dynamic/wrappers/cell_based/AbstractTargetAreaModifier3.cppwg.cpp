#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractTargetAreaModifier.hpp"

#include "AbstractTargetAreaModifier3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTargetAreaModifier<3 > AbstractTargetAreaModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractTargetAreaModifier3_Overrides : public AbstractTargetAreaModifier3{
    public:
    using AbstractTargetAreaModifier3::AbstractTargetAreaModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractTargetAreaModifier3,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractTargetAreaModifier3_class(py::module &m){
py::class_<AbstractTargetAreaModifier3 , AbstractTargetAreaModifier3_Overrides , boost::shared_ptr<AbstractTargetAreaModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "AbstractTargetAreaModifier3")
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractTargetAreaModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &AbstractTargetAreaModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "GetReferenceTargetArea",
            (double(AbstractTargetAreaModifier3::*)()) &AbstractTargetAreaModifier3::GetReferenceTargetArea,
            " "  )
        .def(
            "SetReferenceTargetArea",
            (void(AbstractTargetAreaModifier3::*)(double)) &AbstractTargetAreaModifier3::SetReferenceTargetArea,
            " " , py::arg("referenceTargetArea") )
        .def(
            "UpdateTargetAreas",
            (void(AbstractTargetAreaModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractTargetAreaModifier3::UpdateTargetAreas,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateTargetAreaOfCell",
            (void(AbstractTargetAreaModifier3::*)(::CellPtr const)) &AbstractTargetAreaModifier3::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractTargetAreaModifier3::*)(::out_stream &)) &AbstractTargetAreaModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
