#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "SurfaceAreaConstraintPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef SurfaceAreaConstraintPottsUpdateRule<2 > SurfaceAreaConstraintPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SurfaceAreaConstraintPottsUpdateRule2_Overrides : public SurfaceAreaConstraintPottsUpdateRule2{
    public:
    using SurfaceAreaConstraintPottsUpdateRule2::SurfaceAreaConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            SurfaceAreaConstraintPottsUpdateRule2,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SurfaceAreaConstraintPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_SurfaceAreaConstraintPottsUpdateRule2_class(py::module &m){
py::class_<SurfaceAreaConstraintPottsUpdateRule2 , SurfaceAreaConstraintPottsUpdateRule2_Overrides , boost::shared_ptr<SurfaceAreaConstraintPottsUpdateRule2 >  , AbstractPottsUpdateRule<2>  >(m, "SurfaceAreaConstraintPottsUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(SurfaceAreaConstraintPottsUpdateRule2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &SurfaceAreaConstraintPottsUpdateRule2::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetDeformationEnergyParameter",
            (double(SurfaceAreaConstraintPottsUpdateRule2::*)()) &SurfaceAreaConstraintPottsUpdateRule2::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(SurfaceAreaConstraintPottsUpdateRule2::*)(double)) &SurfaceAreaConstraintPottsUpdateRule2::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "GetMatureCellTargetSurfaceArea",
            (double(SurfaceAreaConstraintPottsUpdateRule2::*)() const ) &SurfaceAreaConstraintPottsUpdateRule2::GetMatureCellTargetSurfaceArea,
            " "  )
        .def(
            "SetMatureCellTargetSurfaceArea",
            (void(SurfaceAreaConstraintPottsUpdateRule2::*)(double)) &SurfaceAreaConstraintPottsUpdateRule2::SetMatureCellTargetSurfaceArea,
            " " , py::arg("matureCellTargetSurfaceArea") )
        .def(
            "OutputUpdateRuleParameters",
            (void(SurfaceAreaConstraintPottsUpdateRule2::*)(::out_stream &)) &SurfaceAreaConstraintPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
