#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"

#include "DifferentialAdhesionGeneralisedLinearSpringForce2_2.cppwg.hpp"

namespace py = pybind11;
typedef DifferentialAdhesionGeneralisedLinearSpringForce<2,2 > DifferentialAdhesionGeneralisedLinearSpringForce2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DifferentialAdhesionGeneralisedLinearSpringForce2_2_Overrides : public DifferentialAdhesionGeneralisedLinearSpringForce2_2{
    public:
    using DifferentialAdhesionGeneralisedLinearSpringForce2_2::DifferentialAdhesionGeneralisedLinearSpringForce;
    double VariableSpringConstantMultiplicationFactor(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<2> & rCellPopulation, bool isCloserThanRestLength) override {
        PYBIND11_OVERRIDE(
            double,
            DifferentialAdhesionGeneralisedLinearSpringForce2_2,
            VariableSpringConstantMultiplicationFactor,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation,
        isCloserThanRestLength);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DifferentialAdhesionGeneralisedLinearSpringForce2_2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_DifferentialAdhesionGeneralisedLinearSpringForce2_2_class(py::module &m){
py::class_<DifferentialAdhesionGeneralisedLinearSpringForce2_2 , DifferentialAdhesionGeneralisedLinearSpringForce2_2_Overrides , boost::shared_ptr<DifferentialAdhesionGeneralisedLinearSpringForce2_2 >  , GeneralisedLinearSpringForce<2>  >(m, "DifferentialAdhesionGeneralisedLinearSpringForce2_2")
        .def(py::init< >())
        .def(
            "VariableSpringConstantMultiplicationFactor",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)(unsigned int, unsigned int, ::AbstractCellPopulation<2> &, bool)) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::VariableSpringConstantMultiplicationFactor,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation"), py::arg("isCloserThanRestLength") )
        .def(
            "GetHomotypicLabelledSpringConstantMultiplier",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)()) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::GetHomotypicLabelledSpringConstantMultiplier,
            " "  )
        .def(
            "SetHomotypicLabelledSpringConstantMultiplier",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)(double)) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::SetHomotypicLabelledSpringConstantMultiplier,
            " " , py::arg("labelledSpringConstantMultiplier") )
        .def(
            "GetHeterotypicSpringConstantMultiplier",
            (double(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)()) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::GetHeterotypicSpringConstantMultiplier,
            " "  )
        .def(
            "SetHeterotypicSpringConstantMultiplier",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)(double)) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::SetHeterotypicSpringConstantMultiplier,
            " " , py::arg("heterotypicSpringConstantMultiplier") )
        .def(
            "OutputForceParameters",
            (void(DifferentialAdhesionGeneralisedLinearSpringForce2_2::*)(::out_stream &)) &DifferentialAdhesionGeneralisedLinearSpringForce2_2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
