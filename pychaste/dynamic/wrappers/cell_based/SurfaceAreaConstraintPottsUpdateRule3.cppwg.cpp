#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "SurfaceAreaConstraintPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef SurfaceAreaConstraintPottsUpdateRule<3 > SurfaceAreaConstraintPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SurfaceAreaConstraintPottsUpdateRule3_Overrides : public SurfaceAreaConstraintPottsUpdateRule3{
    public:
    using SurfaceAreaConstraintPottsUpdateRule3::SurfaceAreaConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            SurfaceAreaConstraintPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SurfaceAreaConstraintPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_SurfaceAreaConstraintPottsUpdateRule3_class(py::module &m){
py::class_<SurfaceAreaConstraintPottsUpdateRule3 , SurfaceAreaConstraintPottsUpdateRule3_Overrides , boost::shared_ptr<SurfaceAreaConstraintPottsUpdateRule3 >  , AbstractPottsUpdateRule<3>  >(m, "SurfaceAreaConstraintPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &SurfaceAreaConstraintPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetDeformationEnergyParameter",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)()) &SurfaceAreaConstraintPottsUpdateRule3::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(double)) &SurfaceAreaConstraintPottsUpdateRule3::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "GetMatureCellTargetSurfaceArea",
            (double(SurfaceAreaConstraintPottsUpdateRule3::*)() const ) &SurfaceAreaConstraintPottsUpdateRule3::GetMatureCellTargetSurfaceArea,
            " "  )
        .def(
            "SetMatureCellTargetSurfaceArea",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(double)) &SurfaceAreaConstraintPottsUpdateRule3::SetMatureCellTargetSurfaceArea,
            " " , py::arg("matureCellTargetSurfaceArea") )
        .def(
            "OutputUpdateRuleParameters",
            (void(SurfaceAreaConstraintPottsUpdateRule3::*)(::out_stream &)) &SurfaceAreaConstraintPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
