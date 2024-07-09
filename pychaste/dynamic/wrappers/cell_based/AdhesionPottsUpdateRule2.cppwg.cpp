#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AdhesionPottsUpdateRule.hpp"

#include "AdhesionPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef AdhesionPottsUpdateRule<2 > AdhesionPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AdhesionPottsUpdateRule2_Overrides : public AdhesionPottsUpdateRule2{
    public:
    using AdhesionPottsUpdateRule2::AdhesionPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule2,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    double GetCellCellAdhesionEnergy(::CellPtr pCellA, ::CellPtr pCellB) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule2,
            GetCellCellAdhesionEnergy,
                    pCellA,
        pCellB);
    }
    double GetCellBoundaryAdhesionEnergy(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            AdhesionPottsUpdateRule2,
            GetCellBoundaryAdhesionEnergy,
                    pCell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AdhesionPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AdhesionPottsUpdateRule2_class(py::module &m){
py::class_<AdhesionPottsUpdateRule2 , AdhesionPottsUpdateRule2_Overrides , boost::shared_ptr<AdhesionPottsUpdateRule2 >  , AbstractPottsUpdateRule<2>  >(m, "AdhesionPottsUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(AdhesionPottsUpdateRule2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &AdhesionPottsUpdateRule2::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "GetCellCellAdhesionEnergy",
            (double(AdhesionPottsUpdateRule2::*)(::CellPtr, ::CellPtr)) &AdhesionPottsUpdateRule2::GetCellCellAdhesionEnergy,
            " " , py::arg("pCellA"), py::arg("pCellB") )
        .def(
            "GetCellBoundaryAdhesionEnergy",
            (double(AdhesionPottsUpdateRule2::*)(::CellPtr)) &AdhesionPottsUpdateRule2::GetCellBoundaryAdhesionEnergy,
            " " , py::arg("pCell") )
        .def(
            "GetCellCellAdhesionEnergyParameter",
            (double(AdhesionPottsUpdateRule2::*)()) &AdhesionPottsUpdateRule2::GetCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetCellBoundaryAdhesionEnergyParameter",
            (double(AdhesionPottsUpdateRule2::*)()) &AdhesionPottsUpdateRule2::GetCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetCellCellAdhesionEnergyParameter",
            (void(AdhesionPottsUpdateRule2::*)(double)) &AdhesionPottsUpdateRule2::SetCellCellAdhesionEnergyParameter,
            " " , py::arg("cellCellAdhesionEnergyEnergyParameter") )
        .def(
            "SetCellBoundaryAdhesionEnergyParameter",
            (void(AdhesionPottsUpdateRule2::*)(double)) &AdhesionPottsUpdateRule2::SetCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("cellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AdhesionPottsUpdateRule2::*)(::out_stream &)) &AdhesionPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
