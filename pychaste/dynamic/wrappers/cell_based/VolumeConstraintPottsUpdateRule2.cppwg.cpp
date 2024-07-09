#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"

#include "VolumeConstraintPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef VolumeConstraintPottsUpdateRule<2 > VolumeConstraintPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeConstraintPottsUpdateRule2_Overrides : public VolumeConstraintPottsUpdateRule2{
    public:
    using VolumeConstraintPottsUpdateRule2::VolumeConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            VolumeConstraintPottsUpdateRule2,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeConstraintPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_VolumeConstraintPottsUpdateRule2_class(py::module &m){
py::class_<VolumeConstraintPottsUpdateRule2 , VolumeConstraintPottsUpdateRule2_Overrides , boost::shared_ptr<VolumeConstraintPottsUpdateRule2 >  , AbstractPottsUpdateRule<2>  >(m, "VolumeConstraintPottsUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(VolumeConstraintPottsUpdateRule2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &VolumeConstraintPottsUpdateRule2::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetDeformationEnergyParameter",
            (double(VolumeConstraintPottsUpdateRule2::*)()) &VolumeConstraintPottsUpdateRule2::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(VolumeConstraintPottsUpdateRule2::*)(double)) &VolumeConstraintPottsUpdateRule2::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "GetMatureCellTargetVolume",
            (double(VolumeConstraintPottsUpdateRule2::*)() const ) &VolumeConstraintPottsUpdateRule2::GetMatureCellTargetVolume,
            " "  )
        .def(
            "SetMatureCellTargetVolume",
            (void(VolumeConstraintPottsUpdateRule2::*)(double)) &VolumeConstraintPottsUpdateRule2::SetMatureCellTargetVolume,
            " " , py::arg("matureCellTargetVolume") )
        .def(
            "OutputUpdateRuleParameters",
            (void(VolumeConstraintPottsUpdateRule2::*)(::out_stream &)) &VolumeConstraintPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
