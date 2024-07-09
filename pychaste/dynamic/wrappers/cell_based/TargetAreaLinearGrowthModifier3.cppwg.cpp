#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"

#include "TargetAreaLinearGrowthModifier3.cppwg.hpp"

namespace py = pybind11;
typedef TargetAreaLinearGrowthModifier<3 > TargetAreaLinearGrowthModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TargetAreaLinearGrowthModifier3_Overrides : public TargetAreaLinearGrowthModifier3{
    public:
    using TargetAreaLinearGrowthModifier3::TargetAreaLinearGrowthModifier;
    void UpdateTargetAreaOfCell(::CellPtr const pCell) override {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier3,
            UpdateTargetAreaOfCell,
                    pCell);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            TargetAreaLinearGrowthModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_TargetAreaLinearGrowthModifier3_class(py::module &m){
py::class_<TargetAreaLinearGrowthModifier3 , TargetAreaLinearGrowthModifier3_Overrides , boost::shared_ptr<TargetAreaLinearGrowthModifier3 >  , AbstractTargetAreaModifier<3>  >(m, "TargetAreaLinearGrowthModifier3")
        .def(py::init< >())
        .def(
            "UpdateTargetAreaOfCell",
            (void(TargetAreaLinearGrowthModifier3::*)(::CellPtr const)) &TargetAreaLinearGrowthModifier3::UpdateTargetAreaOfCell,
            " " , py::arg("pCell") )
        .def(
            "GetAgeToStartGrowing",
            (double(TargetAreaLinearGrowthModifier3::*)()) &TargetAreaLinearGrowthModifier3::GetAgeToStartGrowing,
            " "  )
        .def(
            "SetAgeToStartGrowing",
            (void(TargetAreaLinearGrowthModifier3::*)(double)) &TargetAreaLinearGrowthModifier3::SetAgeToStartGrowing,
            " " , py::arg("ageToStartGrowing") )
        .def(
            "GetGrowthRate",
            (double(TargetAreaLinearGrowthModifier3::*)()) &TargetAreaLinearGrowthModifier3::GetGrowthRate,
            " "  )
        .def(
            "SetGrowthRate",
            (void(TargetAreaLinearGrowthModifier3::*)(double)) &TargetAreaLinearGrowthModifier3::SetGrowthRate,
            " " , py::arg("growthRate") )
        .def(
            "OutputSimulationModifierParameters",
            (void(TargetAreaLinearGrowthModifier3::*)(::out_stream &)) &TargetAreaLinearGrowthModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
