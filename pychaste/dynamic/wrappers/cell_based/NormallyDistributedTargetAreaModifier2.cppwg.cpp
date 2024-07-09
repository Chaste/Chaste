#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NormallyDistributedTargetAreaModifier.hpp"

#include "NormallyDistributedTargetAreaModifier2.cppwg.hpp"

namespace py = pybind11;
typedef NormallyDistributedTargetAreaModifier<2 > NormallyDistributedTargetAreaModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NormallyDistributedTargetAreaModifier2_Overrides : public NormallyDistributedTargetAreaModifier2{
    public:
    using NormallyDistributedTargetAreaModifier2::NormallyDistributedTargetAreaModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            NormallyDistributedTargetAreaModifier2,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NormallyDistributedTargetAreaModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_NormallyDistributedTargetAreaModifier2_class(py::module &m){
py::class_<NormallyDistributedTargetAreaModifier2 , NormallyDistributedTargetAreaModifier2_Overrides , boost::shared_ptr<NormallyDistributedTargetAreaModifier2 >  , AbstractTargetAreaModifier<2>  >(m, "NormallyDistributedTargetAreaModifier2")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(NormallyDistributedTargetAreaModifier2::*)(::CellPtr const)) &NormallyDistributedTargetAreaModifier2::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetGrowthDuration",
            (double(NormallyDistributedTargetAreaModifier2::*)()) &NormallyDistributedTargetAreaModifier2::GetGrowthDuration,
            " "  )
        .def(
            "SetGrowthDuration",
            (void(NormallyDistributedTargetAreaModifier2::*)(double)) &NormallyDistributedTargetAreaModifier2::SetGrowthDuration,
            " " , py::arg("growthDuration") )
        .def(
            "OutputSimulationModifierParameters",
            (void(NormallyDistributedTargetAreaModifier2::*)(::out_stream &)) &NormallyDistributedTargetAreaModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
