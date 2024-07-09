#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractTargetAreaModifier.hpp"

#include "AbstractTargetAreaModifier2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTargetAreaModifier<2 > AbstractTargetAreaModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractTargetAreaModifier2_Overrides : public AbstractTargetAreaModifier2{
    public:
    using AbstractTargetAreaModifier2::AbstractTargetAreaModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractTargetAreaModifier2,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractTargetAreaModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractTargetAreaModifier2_class(py::module &m){
py::class_<AbstractTargetAreaModifier2 , AbstractTargetAreaModifier2_Overrides , boost::shared_ptr<AbstractTargetAreaModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "AbstractTargetAreaModifier2")
        .def(
            "UpdateAtEndOfTimeStep",
            (void(AbstractTargetAreaModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractTargetAreaModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(AbstractTargetAreaModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractTargetAreaModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "GetReferenceTargetArea",
            (double(AbstractTargetAreaModifier2::*)()) &AbstractTargetAreaModifier2::GetReferenceTargetArea,
            " "  )
        .def(
            "SetReferenceTargetArea",
            (void(AbstractTargetAreaModifier2::*)(double)) &AbstractTargetAreaModifier2::SetReferenceTargetArea,
            " " , py::arg("referenceTargetArea") )
        .def(
            "UpdateTargetAreas",
            (void(AbstractTargetAreaModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractTargetAreaModifier2::UpdateTargetAreas,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateTargetAreaOfCell",
            (void(AbstractTargetAreaModifier2::*)(::CellPtr const)) &AbstractTargetAreaModifier2::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractTargetAreaModifier2::*)(::out_stream &)) &AbstractTargetAreaModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
