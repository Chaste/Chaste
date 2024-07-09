#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AdhesionPottsUpdateRule.hpp"

#include "AdhesionPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef AdhesionPottsUpdateRule<3 > AdhesionPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AdhesionPottsUpdateRule3_Overrides : public AdhesionPottsUpdateRule3{
    public:
    using AdhesionPottsUpdateRule3::AdhesionPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    double GetCellCellAdhesionEnergy(::CellPtr pCellA, ::CellPtr pCellB) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule3,
            GetCellCellAdhesionEnergy,
                    pCellA,
        pCellB);
    }
    double GetCellBoundaryAdhesionEnergy(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule3,
            GetCellBoundaryAdhesionEnergy,
                    pCell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AdhesionPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AdhesionPottsUpdateRule3_class(py::module &m){
py::class_<AdhesionPottsUpdateRule3 , AdhesionPottsUpdateRule3_Overrides , boost::shared_ptr<AdhesionPottsUpdateRule3 >  , AbstractPottsUpdateRule<3>  >(m, "AdhesionPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(AdhesionPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &AdhesionPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetCellCellAdhesionEnergy",
            (double(AdhesionPottsUpdateRule3::*)(::CellPtr, ::CellPtr)) &AdhesionPottsUpdateRule3::GetCellCellAdhesionEnergy,
            " " , py::arg("pCellA"), py::arg("pCellB") )
        .def(
            "GetCellBoundaryAdhesionEnergy",
            (double(AdhesionPottsUpdateRule3::*)(::CellPtr)) &AdhesionPottsUpdateRule3::GetCellBoundaryAdhesionEnergy,
            " " , py::arg("pCell") )
        .def(
            "GetCellCellAdhesionEnergyParameter",
            (double(AdhesionPottsUpdateRule3::*)()) &AdhesionPottsUpdateRule3::GetCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetCellBoundaryAdhesionEnergyParameter",
            (double(AdhesionPottsUpdateRule3::*)()) &AdhesionPottsUpdateRule3::GetCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetCellCellAdhesionEnergyParameter",
            (void(AdhesionPottsUpdateRule3::*)(double)) &AdhesionPottsUpdateRule3::SetCellCellAdhesionEnergyParameter,
            " " , py::arg("cellCellAdhesionEnergyEnergyParameter") )
        .def(
            "SetCellBoundaryAdhesionEnergyParameter",
            (void(AdhesionPottsUpdateRule3::*)(double)) &AdhesionPottsUpdateRule3::SetCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("cellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AdhesionPottsUpdateRule3::*)(::out_stream &)) &AdhesionPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
