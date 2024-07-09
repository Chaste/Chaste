#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "SimpleTargetAreaModifier3.cppwg.hpp"

namespace py = pybind11;
typedef SimpleTargetAreaModifier<3 > SimpleTargetAreaModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SimpleTargetAreaModifier3_Overrides : public SimpleTargetAreaModifier3{
    public:
    using SimpleTargetAreaModifier3::SimpleTargetAreaModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            SimpleTargetAreaModifier3,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SimpleTargetAreaModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_SimpleTargetAreaModifier3_class(py::module &m){
py::class_<SimpleTargetAreaModifier3 , SimpleTargetAreaModifier3_Overrides , boost::shared_ptr<SimpleTargetAreaModifier3 >  , AbstractTargetAreaModifier<3>  >(m, "SimpleTargetAreaModifier3")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(SimpleTargetAreaModifier3::*)(::CellPtr const)) &SimpleTargetAreaModifier3::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetGrowthDuration",
            (double(SimpleTargetAreaModifier3::*)()) &SimpleTargetAreaModifier3::GetGrowthDuration,
            " "  )
        .def(
            "SetGrowthDuration",
            (void(SimpleTargetAreaModifier3::*)(double)) &SimpleTargetAreaModifier3::SetGrowthDuration,
            " " , py::arg("growthDuration") )
        .def(
            "OutputSimulationModifierParameters",
            (void(SimpleTargetAreaModifier3::*)(::out_stream &)) &SimpleTargetAreaModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
