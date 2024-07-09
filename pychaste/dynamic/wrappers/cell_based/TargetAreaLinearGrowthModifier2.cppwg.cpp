#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"

#include "TargetAreaLinearGrowthModifier2.cppwg.hpp"

namespace py = pybind11;
typedef TargetAreaLinearGrowthModifier<2 > TargetAreaLinearGrowthModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TargetAreaLinearGrowthModifier2_Overrides : public TargetAreaLinearGrowthModifier2{
    public:
    using TargetAreaLinearGrowthModifier2::TargetAreaLinearGrowthModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier2,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_TargetAreaLinearGrowthModifier2_class(py::module &m){
py::class_<TargetAreaLinearGrowthModifier2 , TargetAreaLinearGrowthModifier2_Overrides , boost::shared_ptr<TargetAreaLinearGrowthModifier2 >  , AbstractTargetAreaModifier<2>  >(m, "TargetAreaLinearGrowthModifier2")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(TargetAreaLinearGrowthModifier2::*)(::CellPtr const)) &TargetAreaLinearGrowthModifier2::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetAgeToStartGrowing",
            (double(TargetAreaLinearGrowthModifier2::*)()) &TargetAreaLinearGrowthModifier2::GetAgeToStartGrowing,
            " "  )
        .def(
            "SetAgeToStartGrowing",
            (void(TargetAreaLinearGrowthModifier2::*)(double)) &TargetAreaLinearGrowthModifier2::SetAgeToStartGrowing,
            " " , py::arg("ageToStartGrowing") )
        .def(
            "GetGrowthRate",
            (double(TargetAreaLinearGrowthModifier2::*)()) &TargetAreaLinearGrowthModifier2::GetGrowthRate,
            " "  )
        .def(
            "SetGrowthRate",
            (void(TargetAreaLinearGrowthModifier2::*)(double)) &TargetAreaLinearGrowthModifier2::SetGrowthRate,
            " " , py::arg("growthRate") )
        .def(
            "OutputSimulationModifierParameters",
            (void(TargetAreaLinearGrowthModifier2::*)(::out_stream &)) &TargetAreaLinearGrowthModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
