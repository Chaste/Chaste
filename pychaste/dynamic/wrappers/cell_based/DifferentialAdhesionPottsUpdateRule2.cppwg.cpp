#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "DifferentialAdhesionPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef DifferentialAdhesionPottsUpdateRule<2 > DifferentialAdhesionPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DifferentialAdhesionPottsUpdateRule2_Overrides : public DifferentialAdhesionPottsUpdateRule2{
    public:
    using DifferentialAdhesionPottsUpdateRule2::DifferentialAdhesionPottsUpdateRule;
    double GetCellCellAdhesionEnergy(::CellPtr pCellA, ::CellPtr pCellB) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule2,
            GetCellCellAdhesionEnergy,
                    pCellA,
        pCellB);
    }
    double GetCellBoundaryAdhesionEnergy(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule2,
            GetCellBoundaryAdhesionEnergy,
                    pCell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DifferentialAdhesionPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_DifferentialAdhesionPottsUpdateRule2_class(py::module &m){
py::class_<DifferentialAdhesionPottsUpdateRule2 , DifferentialAdhesionPottsUpdateRule2_Overrides , boost::shared_ptr<DifferentialAdhesionPottsUpdateRule2 >  , AdhesionPottsUpdateRule<2>  >(m, "DifferentialAdhesionPottsUpdateRule2")
        .def(py::init< >())
        .def(
            "GetCellCellAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule2::*)(::CellPtr, ::CellPtr)) &DifferentialAdhesionPottsUpdateRule2::GetCellCellAdhesionEnergy,
            " " , py::arg("pCellA"), py::arg("pCellB") )
        .def(
            "GetCellBoundaryAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule2::*)(::CellPtr)) &DifferentialAdhesionPottsUpdateRule2::GetCellBoundaryAdhesionEnergy,
            " " , py::arg("pCell") )
        .def(
            "GetLabelledCellLabelledCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule2::*)()) &DifferentialAdhesionPottsUpdateRule2::GetLabelledCellLabelledCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetLabelledCellCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule2::*)()) &DifferentialAdhesionPottsUpdateRule2::GetLabelledCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetLabelledCellBoundaryAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule2::*)()) &DifferentialAdhesionPottsUpdateRule2::GetLabelledCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetLabelledCellLabelledCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule2::*)(double)) &DifferentialAdhesionPottsUpdateRule2::SetLabelledCellLabelledCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellLabelledCellAdhesionEnergyParameter") )
        .def(
            "SetLabelledCellCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule2::*)(double)) &DifferentialAdhesionPottsUpdateRule2::SetLabelledCellCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellCellAdhesionEnergyParameter") )
        .def(
            "SetLabelledCellBoundaryAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule2::*)(double)) &DifferentialAdhesionPottsUpdateRule2::SetLabelledCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("labelledCellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(DifferentialAdhesionPottsUpdateRule2::*)(::out_stream &)) &DifferentialAdhesionPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
