#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"

#include "DifferentialAdhesionGeneralisedLinearSpringForce3_3.cppwg.hpp"

namespace py = pybind11;
typedef DifferentialAdhesionGeneralisedLinearSpringForce<3,3 > DifferentialAdhesionGeneralisedLinearSpringForce3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DifferentialAdhesionGeneralisedLinearSpringForce3_3_Overrides : public DifferentialAdhesionGeneralisedLinearSpringForce3_3{
    public:
    using DifferentialAdhesionGeneralisedLinearSpringForce3_3::DifferentialAdhesionGeneralisedLinearSpringForce;
    double VariableSpringConstantMultiplicationFactor(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation, bool isCloserThanRestLength) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionGeneralisedLinearSpringForce3_3,
            VariableSpringConstantMultiplicationFactor,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation,
        isCloserThanRestLength);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DifferentialAdhesionGeneralisedLinearSpringForce3_3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_DifferentialAdhesionGeneralisedLinearSpringForce3_3_class(py::module &m){
py::class_<DifferentialAdhesionGeneralisedLinearSpringForce3_3 , DifferentialAdhesionGeneralisedLinearSpringForce3_3_Overrides , boost::shared_ptr<DifferentialAdhesionGeneralisedLinearSpringForce3_3 >  , GeneralisedLinearSpringForce<3>  >(m, "DifferentialAdhesionGeneralisedLinearSpringForce3_3")
        .def(py::init< >())
        .def(
            "VariableSpringConstantMultiplicationFactor",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &, bool)) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::VariableSpringConstantMultiplicationFactor,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation"), py::arg("isCloserThanRestLength") )
        .def(
            "GetHomotypicLabelledSpringConstantMultiplier",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)()) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::GetHomotypicLabelledSpringConstantMultiplier,
            " "  )
        .def(
            "SetHomotypicLabelledSpringConstantMultiplier",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)(double)) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::SetHomotypicLabelledSpringConstantMultiplier,
            " " , py::arg("labelledSpringConstantMultiplier") )
        .def(
            "GetHeterotypicSpringConstantMultiplier",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)()) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::GetHeterotypicSpringConstantMultiplier,
            " "  )
        .def(
            "SetHeterotypicSpringConstantMultiplier",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)(double)) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::SetHeterotypicSpringConstantMultiplier,
            " " , py::arg("heterotypicSpringConstantMultiplier") )
        .def(
            "OutputForceParameters",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce3_3::*)(::out_stream &)) &DifferentialAdhesionGeneralisedLinearSpringForce3_3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
