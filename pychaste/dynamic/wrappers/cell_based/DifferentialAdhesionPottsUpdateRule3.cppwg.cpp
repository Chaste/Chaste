#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "DifferentialAdhesionPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef DifferentialAdhesionPottsUpdateRule<3 > DifferentialAdhesionPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DifferentialAdhesionPottsUpdateRule3_Overrides : public DifferentialAdhesionPottsUpdateRule3{
    public:
    using DifferentialAdhesionPottsUpdateRule3::DifferentialAdhesionPottsUpdateRule;
    double GetCellCellAdhesionEnergy(::CellPtr pCellA, ::CellPtr pCellB) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule3,
            GetCellCellAdhesionEnergy,
                    pCellA,
        pCellB);
    }
    double GetCellBoundaryAdhesionEnergy(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionPottsUpdateRule3,
            GetCellBoundaryAdhesionEnergy,
                    pCell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DifferentialAdhesionPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_DifferentialAdhesionPottsUpdateRule3_class(py::module &m){
py::class_<DifferentialAdhesionPottsUpdateRule3 , DifferentialAdhesionPottsUpdateRule3_Overrides , boost::shared_ptr<DifferentialAdhesionPottsUpdateRule3 >  , AdhesionPottsUpdateRule<3>  >(m, "DifferentialAdhesionPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "GetCellCellAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule3::*)(::CellPtr, ::CellPtr)) &DifferentialAdhesionPottsUpdateRule3::GetCellCellAdhesionEnergy,
            " " , py::arg("pCellA"), py::arg("pCellB") )
        .def(
            "GetCellBoundaryAdhesionEnergy",
            (double(DifferentialAdhesionPottsUpdateRule3::*)(::CellPtr)) &DifferentialAdhesionPottsUpdateRule3::GetCellBoundaryAdhesionEnergy,
            " " , py::arg("pCell") )
        .def(
            "GetLabelledCellLabelledCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule3::*)()) &DifferentialAdhesionPottsUpdateRule3::GetLabelledCellLabelledCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetLabelledCellCellAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule3::*)()) &DifferentialAdhesionPottsUpdateRule3::GetLabelledCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetLabelledCellBoundaryAdhesionEnergyParameter",
            (double(DifferentialAdhesionPottsUpdateRule3::*)()) &DifferentialAdhesionPottsUpdateRule3::GetLabelledCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetLabelledCellLabelledCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule3::*)(double)) &DifferentialAdhesionPottsUpdateRule3::SetLabelledCellLabelledCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellLabelledCellAdhesionEnergyParameter") )
        .def(
            "SetLabelledCellCellAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule3::*)(double)) &DifferentialAdhesionPottsUpdateRule3::SetLabelledCellCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellCellAdhesionEnergyParameter") )
        .def(
            "SetLabelledCellBoundaryAdhesionEnergyParameter",
            (void(DifferentialAdhesionPottsUpdateRule3::*)(double)) &DifferentialAdhesionPottsUpdateRule3::SetLabelledCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("labelledCellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(DifferentialAdhesionPottsUpdateRule3::*)(::out_stream &)) &DifferentialAdhesionPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
