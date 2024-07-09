#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include "NagaiHondaDifferentialAdhesionForce3.cppwg.hpp"

namespace py = pybind11;
typedef NagaiHondaDifferentialAdhesionForce<3 > NagaiHondaDifferentialAdhesionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NagaiHondaDifferentialAdhesionForce3_Overrides : public NagaiHondaDifferentialAdhesionForce3{
    public:
    using NagaiHondaDifferentialAdhesionForce3::NagaiHondaDifferentialAdhesionForce;
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaDifferentialAdhesionForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_NagaiHondaDifferentialAdhesionForce3_class(py::module &m){
py::class_<NagaiHondaDifferentialAdhesionForce3 , NagaiHondaDifferentialAdhesionForce3_Overrides , boost::shared_ptr<NagaiHondaDifferentialAdhesionForce3 >  , NagaiHondaForce<3>  >(m, "NagaiHondaDifferentialAdhesionForce3")
        .def(py::init< >())
        .def(
            "GetNagaiHondaLabelledCellCellAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce3::*)()) &NagaiHondaDifferentialAdhesionForce3::GetNagaiHondaLabelledCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce3::*)()) &NagaiHondaDifferentialAdhesionForce3::GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce3::*)()) &NagaiHondaDifferentialAdhesionForce3::GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetNagaiHondaLabelledCellCellAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce3::*)(double)) &NagaiHondaDifferentialAdhesionForce3::SetNagaiHondaLabelledCellCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellCellAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce3::*)(double)) &NagaiHondaDifferentialAdhesionForce3::SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellLabelledCellAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce3::*)(double)) &NagaiHondaDifferentialAdhesionForce3::SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("labelledCellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputForceParameters",
            (void(NagaiHondaDifferentialAdhesionForce3::*)(::out_stream &)) &NagaiHondaDifferentialAdhesionForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
