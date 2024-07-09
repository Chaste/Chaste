#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FarhadifarForce.hpp"

#include "FarhadifarForce2.cppwg.hpp"

namespace py = pybind11;
typedef FarhadifarForce<2 > FarhadifarForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class FarhadifarForce2_Overrides : public FarhadifarForce2{
    public:
    using FarhadifarForce2::FarhadifarForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    double GetLineTensionParameter(::Node<2> * pNodeA, ::Node<2> * pNodeB, ::VertexBasedCellPopulation<2> & rVertexCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            FarhadifarForce2,
            GetLineTensionParameter,
                    pNodeA,
        pNodeB,
        rVertexCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            FarhadifarForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_FarhadifarForce2_class(py::module &m){
py::class_<FarhadifarForce2 , FarhadifarForce2_Overrides , boost::shared_ptr<FarhadifarForce2 >  , AbstractForce<2>  >(m, "FarhadifarForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(FarhadifarForce2::*)(::AbstractCellPopulation<2> &)) &FarhadifarForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetLineTensionParameter",
            (double(FarhadifarForce2::*)(::Node<2> *, ::Node<2> *, ::VertexBasedCellPopulation<2> &)) &FarhadifarForce2::GetLineTensionParameter,
            " " , py::arg("pNodeA"), py::arg("pNodeB"), py::arg("rVertexCellPopulation") )
        .def(
            "GetAreaElasticityParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetAreaElasticityParameter,
            " "  )
        .def(
            "GetPerimeterContractilityParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetPerimeterContractilityParameter,
            " "  )
        .def(
            "GetLineTensionParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetLineTensionParameter,
            " "  )
        .def(
            "GetBoundaryLineTensionParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetBoundaryLineTensionParameter,
            " "  )
        .def(
            "GetTargetAreaParameter",
            (double(FarhadifarForce2::*)()) &FarhadifarForce2::GetTargetAreaParameter,
            " "  )
        .def(
            "SetAreaElasticityParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetAreaElasticityParameter,
            " " , py::arg("areaElasticityParameter") )
        .def(
            "SetPerimeterContractilityParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetPerimeterContractilityParameter,
            " " , py::arg("perimeterContractilityParameter") )
        .def(
            "SetLineTensionParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetLineTensionParameter,
            " " , py::arg("lineTensionParameter") )
        .def(
            "SetBoundaryLineTensionParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetBoundaryLineTensionParameter,
            " " , py::arg("boundaryLineTensionParameter") )
        .def(
            "SetTargetAreaParameter",
            (void(FarhadifarForce2::*)(double)) &FarhadifarForce2::SetTargetAreaParameter,
            " " , py::arg("targetAreaParameter") )
        .def(
            "OutputForceParameters",
            (void(FarhadifarForce2::*)(::out_stream &)) &FarhadifarForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
