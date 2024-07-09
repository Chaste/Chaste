#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ExtrinsicPullModifier.hpp"

#include "ExtrinsicPullModifier2.cppwg.hpp"

namespace py = pybind11;
typedef ExtrinsicPullModifier<2 > ExtrinsicPullModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ExtrinsicPullModifier2_Overrides : public ExtrinsicPullModifier2{
    public:
    using ExtrinsicPullModifier2::ExtrinsicPullModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ExtrinsicPullModifier2_class(py::module &m){
py::class_<ExtrinsicPullModifier2 , ExtrinsicPullModifier2_Overrides , boost::shared_ptr<ExtrinsicPullModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "ExtrinsicPullModifier2")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ExtrinsicPullModifier2::*)(::AbstractCellPopulation<2> &)) &ExtrinsicPullModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ExtrinsicPullModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ExtrinsicPullModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetApplyExtrinsicPullToAllNodes",
            (void(ExtrinsicPullModifier2::*)(bool)) &ExtrinsicPullModifier2::SetApplyExtrinsicPullToAllNodes,
            " " , py::arg("applyExtrinsicPullToAllNodes") )
        .def(
            "SetSpeed",
            (void(ExtrinsicPullModifier2::*)(double)) &ExtrinsicPullModifier2::SetSpeed,
            " " , py::arg("speed") )
        .def(
            "GetApplyExtrinsicPullToAllNodes",
            (bool(ExtrinsicPullModifier2::*)()) &ExtrinsicPullModifier2::GetApplyExtrinsicPullToAllNodes,
            " "  )
        .def(
            "GetSpeed",
            (double(ExtrinsicPullModifier2::*)()) &ExtrinsicPullModifier2::GetSpeed,
            " "  )
        .def(
            "OutputSimulationModifierParameters",
            (void(ExtrinsicPullModifier2::*)(::out_stream &)) &ExtrinsicPullModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
