#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "GeneralisedLinearSpringForce3_3.cppwg.hpp"

namespace py = pybind11;
typedef GeneralisedLinearSpringForce<3,3 > GeneralisedLinearSpringForce3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class GeneralisedLinearSpringForce3_3_Overrides : public GeneralisedLinearSpringForce3_3{
    public:
    using GeneralisedLinearSpringForce3_3::GeneralisedLinearSpringForce;
    double VariableSpringConstantMultiplicationFactor(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation, bool isCloserThanRestLength) override {
        PYBIND11_OVERRIDE(
            double,
            GeneralisedLinearSpringForce3_3,
            VariableSpringConstantMultiplicationFactor,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation,
        isCloserThanRestLength);
    }
    ::boost::numeric::ublas::c_vector<double, 3> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            GeneralisedLinearSpringForce3_3,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            GeneralisedLinearSpringForce3_3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_GeneralisedLinearSpringForce3_3_class(py::module &m){
py::class_<GeneralisedLinearSpringForce3_3 , GeneralisedLinearSpringForce3_3_Overrides , boost::shared_ptr<GeneralisedLinearSpringForce3_3 >  , AbstractTwoBodyInteractionForce<3>  >(m, "GeneralisedLinearSpringForce3_3")
        .def(py::init< >())
        .def(
            "VariableSpringConstantMultiplicationFactor",
            (double(GeneralisedLinearSpringForce3_3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &, bool)) &GeneralisedLinearSpringForce3_3::VariableSpringConstantMultiplicationFactor,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation"), py::arg("isCloserThanRestLength") )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 3>(GeneralisedLinearSpringForce3_3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &)) &GeneralisedLinearSpringForce3_3::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "GetMeinekeSpringStiffness",
            (double(GeneralisedLinearSpringForce3_3::*)()) &GeneralisedLinearSpringForce3_3::GetMeinekeSpringStiffness,
            " "  )
        .def(
            "GetMeinekeDivisionRestingSpringLength",
            (double(GeneralisedLinearSpringForce3_3::*)()) &GeneralisedLinearSpringForce3_3::GetMeinekeDivisionRestingSpringLength,
            " "  )
        .def(
            "GetMeinekeSpringGrowthDuration",
            (double(GeneralisedLinearSpringForce3_3::*)()) &GeneralisedLinearSpringForce3_3::GetMeinekeSpringGrowthDuration,
            " "  )
        .def(
            "SetMeinekeSpringStiffness",
            (void(GeneralisedLinearSpringForce3_3::*)(double)) &GeneralisedLinearSpringForce3_3::SetMeinekeSpringStiffness,
            " " , py::arg("springStiffness") )
        .def(
            "SetMeinekeDivisionRestingSpringLength",
            (void(GeneralisedLinearSpringForce3_3::*)(double)) &GeneralisedLinearSpringForce3_3::SetMeinekeDivisionRestingSpringLength,
            " " , py::arg("divisionRestingSpringLength") )
        .def(
            "SetMeinekeSpringGrowthDuration",
            (void(GeneralisedLinearSpringForce3_3::*)(double)) &GeneralisedLinearSpringForce3_3::SetMeinekeSpringGrowthDuration,
            " " , py::arg("springGrowthDuration") )
        .def(
            "OutputForceParameters",
            (void(GeneralisedLinearSpringForce3_3::*)(::out_stream &)) &GeneralisedLinearSpringForce3_3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
