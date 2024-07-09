#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NagaiHondaForce.hpp"

#include "NagaiHondaForce3.cppwg.hpp"

namespace py = pybind11;
typedef NagaiHondaForce<3 > NagaiHondaForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NagaiHondaForce3_Overrides : public NagaiHondaForce3{
    public:
    using NagaiHondaForce3::NagaiHondaForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    double GetAdhesionParameter(::Node<3> * pNodeA, ::Node<3> * pNodeB, ::VertexBasedCellPopulation<3> & rVertexCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            NagaiHondaForce3,
            GetAdhesionParameter,
                    pNodeA,
        pNodeB,
        rVertexCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NagaiHondaForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_NagaiHondaForce3_class(py::module &m){
py::class_<NagaiHondaForce3 , NagaiHondaForce3_Overrides , boost::shared_ptr<NagaiHondaForce3 >  , AbstractForce<3>  >(m, "NagaiHondaForce3")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(NagaiHondaForce3::*)(::AbstractCellPopulation<3> &)) &NagaiHondaForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetAdhesionParameter",
            (double(NagaiHondaForce3::*)(::Node<3> *, ::Node<3> *, ::VertexBasedCellPopulation<3> &)) &NagaiHondaForce3::GetAdhesionParameter,
            " " , py::arg("pNodeA"), py::arg("pNodeB"), py::arg("rVertexCellPopulation") )
        .def(
            "GetNagaiHondaDeformationEnergyParameter",
            (double(NagaiHondaForce3::*)()) &NagaiHondaForce3::GetNagaiHondaDeformationEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaMembraneSurfaceEnergyParameter",
            (double(NagaiHondaForce3::*)()) &NagaiHondaForce3::GetNagaiHondaMembraneSurfaceEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaCellCellAdhesionEnergyParameter",
            (double(NagaiHondaForce3::*)()) &NagaiHondaForce3::GetNagaiHondaCellCellAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaCellBoundaryAdhesionEnergyParameter",
            (double(NagaiHondaForce3::*)()) &NagaiHondaForce3::GetNagaiHondaCellBoundaryAdhesionEnergyParameter,
            " "  )
        .def(
            "GetNagaiHondaTargetAreaParameter",
            (double(NagaiHondaForce3::*)()) &NagaiHondaForce3::GetNagaiHondaTargetAreaParameter,
            " "  )
        .def(
            "SetNagaiHondaDeformationEnergyParameter",
            (void(NagaiHondaForce3::*)(double)) &NagaiHondaForce3::SetNagaiHondaDeformationEnergyParameter,
            " " , py::arg("nagaiHondaDeformationEnergyParameter") )
        .def(
            "SetNagaiHondaMembraneSurfaceEnergyParameter",
            (void(NagaiHondaForce3::*)(double)) &NagaiHondaForce3::SetNagaiHondaMembraneSurfaceEnergyParameter,
            " " , py::arg("nagaiHondaMembraneSurfaceEnergyParameter") )
        .def(
            "SetNagaiHondaCellCellAdhesionEnergyParameter",
            (void(NagaiHondaForce3::*)(double)) &NagaiHondaForce3::SetNagaiHondaCellCellAdhesionEnergyParameter,
            " " , py::arg("nagaiHondaCellCellAdhesionEnergyEnergyParameter") )
        .def(
            "SetNagaiHondaCellBoundaryAdhesionEnergyParameter",
            (void(NagaiHondaForce3::*)(double)) &NagaiHondaForce3::SetNagaiHondaCellBoundaryAdhesionEnergyParameter,
            " " , py::arg("nagaiHondaCellBoundaryAdhesionEnergyParameter") )
        .def(
            "SetNagaiHondaTargetAreaParameter",
            (void(NagaiHondaForce3::*)(double)) &NagaiHondaForce3::SetNagaiHondaTargetAreaParameter,
            " " , py::arg("nagaiHondaTargetAreaParameter") )
        .def(
            "OutputForceParameters",
            (void(NagaiHondaForce3::*)(::out_stream &)) &NagaiHondaForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
