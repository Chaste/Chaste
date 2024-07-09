#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ExtrinsicPullModifier.hpp"

#include "ExtrinsicPullModifier3.cppwg.hpp"

namespace py = pybind11;
typedef ExtrinsicPullModifier<3 > ExtrinsicPullModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ExtrinsicPullModifier3_Overrides : public ExtrinsicPullModifier3{
    public:
    using ExtrinsicPullModifier3::ExtrinsicPullModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ExtrinsicPullModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ExtrinsicPullModifier3_class(py::module &m){
py::class_<ExtrinsicPullModifier3 , ExtrinsicPullModifier3_Overrides , boost::shared_ptr<ExtrinsicPullModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "ExtrinsicPullModifier3")
        .def(py::init< >())
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ExtrinsicPullModifier3::*)(::AbstractCellPopulation<3> &)) &ExtrinsicPullModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ExtrinsicPullModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &ExtrinsicPullModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetApplyExtrinsicPullToAllNodes",
            (void(ExtrinsicPullModifier3::*)(bool)) &ExtrinsicPullModifier3::SetApplyExtrinsicPullToAllNodes,
            " " , py::arg("applyExtrinsicPullToAllNodes") )
        .def(
            "SetSpeed",
            (void(ExtrinsicPullModifier3::*)(double)) &ExtrinsicPullModifier3::SetSpeed,
            " " , py::arg("speed") )
        .def(
            "GetApplyExtrinsicPullToAllNodes",
            (bool(ExtrinsicPullModifier3::*)()) &ExtrinsicPullModifier3::GetApplyExtrinsicPullToAllNodes,
            " "  )
        .def(
            "GetSpeed",
            (double(ExtrinsicPullModifier3::*)()) &ExtrinsicPullModifier3::GetSpeed,
            " "  )
        .def(
            "OutputSimulationModifierParameters",
            (void(ExtrinsicPullModifier3::*)(::out_stream &)) &ExtrinsicPullModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
