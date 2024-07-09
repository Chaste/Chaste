#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include "NagaiHondaDifferentialAdhesionForce2.cppwg.hpp"

namespace py = pybind11;
typedef NagaiHondaDifferentialAdhesionForce<2 > NagaiHondaDifferentialAdhesionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NagaiHondaDifferentialAdhesionForce2_Overrides : public NagaiHondaDifferentialAdhesionForce2{
    public:
    using NagaiHondaDifferentialAdhesionForce2::NagaiHondaDifferentialAdhesionForce;
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaDifferentialAdhesionForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_NagaiHondaDifferentialAdhesionForce2_class(py::module &m){
py::class_<NagaiHondaDifferentialAdhesionForce2 , NagaiHondaDifferentialAdhesionForce2_Overrides , boost::shared_ptr<NagaiHondaDifferentialAdhesionForce2 >  , NagaiHondaForce<2>  >(m, "NagaiHondaDifferentialAdhesionForce2")
        .def(py::init< >())
        .def(
            "GetNagaiHondaLabelledCellCellAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce2::*)()) &NagaiHondaDifferentialAdhesionForce2::GetNagaiHondaLabelledCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce2::*)()) &NagaiHondaDifferentialAdhesionForce2::GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter",
            (double(NagaiHondaDifferentialAdhesionForce2::*)()) &NagaiHondaDifferentialAdhesionForce2::GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "SetNagaiHondaLabelledCellCellAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce2::*)(double)) &NagaiHondaDifferentialAdhesionForce2::SetNagaiHondaLabelledCellCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellCellAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce2::*)(double)) &NagaiHondaDifferentialAdhesionForce2::SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter,
            " " , py::arg("labelledCellLabelledCellAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter",
            (void(NagaiHondaDifferentialAdhesionForce2::*)(double)) &NagaiHondaDifferentialAdhesionForce2::SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("labelledCellBoundaryAdhesionEnergyParameter") )
        .def(
            "OutputForceParameters",
            (void(NagaiHondaDifferentialAdhesionForce2::*)(::out_stream &)) &NagaiHondaDifferentialAdhesionForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
