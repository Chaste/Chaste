#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NormallyDistributedTargetAreaModifier.hpp"

#include "NormallyDistributedTargetAreaModifier3.cppwg.hpp"

namespace py = pybind11;
typedef NormallyDistributedTargetAreaModifier<3 > NormallyDistributedTargetAreaModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NormallyDistributedTargetAreaModifier3_Overrides : public NormallyDistributedTargetAreaModifier3{
    public:
    using NormallyDistributedTargetAreaModifier3::NormallyDistributedTargetAreaModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            NormallyDistributedTargetAreaModifier3,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NormallyDistributedTargetAreaModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_NormallyDistributedTargetAreaModifier3_class(py::module &m){
py::class_<NormallyDistributedTargetAreaModifier3 , NormallyDistributedTargetAreaModifier3_Overrides , boost::shared_ptr<NormallyDistributedTargetAreaModifier3 >  , AbstractTargetAreaModifier<3>  >(m, "NormallyDistributedTargetAreaModifier3")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(NormallyDistributedTargetAreaModifier3::*)(::CellPtr const)) &NormallyDistributedTargetAreaModifier3::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetGrowthDuration",
            (double(NormallyDistributedTargetAreaModifier3::*)()) &NormallyDistributedTargetAreaModifier3::GetGrowthDuration,
            " "  )
        .def(
            "SetGrowthDuration",
            (void(NormallyDistributedTargetAreaModifier3::*)(double)) &NormallyDistributedTargetAreaModifier3::SetGrowthDuration,
            " " , py::arg("growthDuration") )
        .def(
            "OutputSimulationModifierParameters",
            (void(NormallyDistributedTargetAreaModifier3::*)(::out_stream &)) &NormallyDistributedTargetAreaModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
