#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "SimpleTargetAreaModifier2.cppwg.hpp"

namespace py = pybind11;
typedef SimpleTargetAreaModifier<2 > SimpleTargetAreaModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SimpleTargetAreaModifier2_Overrides : public SimpleTargetAreaModifier2{
    public:
    using SimpleTargetAreaModifier2::SimpleTargetAreaModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            SimpleTargetAreaModifier2,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SimpleTargetAreaModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_SimpleTargetAreaModifier2_class(py::module &m){
py::class_<SimpleTargetAreaModifier2 , SimpleTargetAreaModifier2_Overrides , boost::shared_ptr<SimpleTargetAreaModifier2 >  , AbstractTargetAreaModifier<2>  >(m, "SimpleTargetAreaModifier2")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(SimpleTargetAreaModifier2::*)(::CellPtr const)) &SimpleTargetAreaModifier2::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetGrowthDuration",
            (double(SimpleTargetAreaModifier2::*)()) &SimpleTargetAreaModifier2::GetGrowthDuration,
            " "  )
        .def(
            "SetGrowthDuration",
            (void(SimpleTargetAreaModifier2::*)(double)) &SimpleTargetAreaModifier2::SetGrowthDuration,
            " " , py::arg("growthDuration") )
        .def(
            "OutputSimulationModifierParameters",
            (void(SimpleTargetAreaModifier2::*)(::out_stream &)) &SimpleTargetAreaModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
