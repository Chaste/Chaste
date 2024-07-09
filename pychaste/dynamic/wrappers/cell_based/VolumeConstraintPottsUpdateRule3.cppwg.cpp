#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"

#include "VolumeConstraintPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef VolumeConstraintPottsUpdateRule<3 > VolumeConstraintPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeConstraintPottsUpdateRule3_Overrides : public VolumeConstraintPottsUpdateRule3{
    public:
    using VolumeConstraintPottsUpdateRule3::VolumeConstraintPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            VolumeConstraintPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeConstraintPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_VolumeConstraintPottsUpdateRule3_class(py::module &m){
py::class_<VolumeConstraintPottsUpdateRule3 , VolumeConstraintPottsUpdateRule3_Overrides , boost::shared_ptr<VolumeConstraintPottsUpdateRule3 >  , AbstractPottsUpdateRule<3>  >(m, "VolumeConstraintPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(VolumeConstraintPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &VolumeConstraintPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetDeformationEnergyParameter",
            (double(VolumeConstraintPottsUpdateRule3::*)()) &VolumeConstraintPottsUpdateRule3::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(VolumeConstraintPottsUpdateRule3::*)(double)) &VolumeConstraintPottsUpdateRule3::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "GetMatureCellTargetVolume",
            (double(VolumeConstraintPottsUpdateRule3::*)() const ) &VolumeConstraintPottsUpdateRule3::GetMatureCellTargetVolume,
            " "  )
        .def(
            "SetMatureCellTargetVolume",
            (void(VolumeConstraintPottsUpdateRule3::*)(double)) &VolumeConstraintPottsUpdateRule3::SetMatureCellTargetVolume,
            " " , py::arg("matureCellTargetVolume") )
        .def(
            "OutputUpdateRuleParameters",
            (void(VolumeConstraintPottsUpdateRule3::*)(::out_stream &)) &VolumeConstraintPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
