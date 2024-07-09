#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DiffusionCaUpdateRule.hpp"

#include "DiffusionCaUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef DiffusionCaUpdateRule<2 > DiffusionCaUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DiffusionCaUpdateRule2_Overrides : public DiffusionCaUpdateRule2{
    public:
    using DiffusionCaUpdateRule2::DiffusionCaUpdateRule;
    double EvaluateProbability(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::CaBasedCellPopulation<2> & rCellPopulation, double dt, double deltaX, ::CellPtr cell) override {
        PYBIND11_OVERRIDE(
            double,
            DiffusionCaUpdateRule2,
            EvaluateProbability,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation,
        dt,
        deltaX,
        cell);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DiffusionCaUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_DiffusionCaUpdateRule2_class(py::module &m){
py::class_<DiffusionCaUpdateRule2 , DiffusionCaUpdateRule2_Overrides , boost::shared_ptr<DiffusionCaUpdateRule2 >  , AbstractCaUpdateRule<2>  >(m, "DiffusionCaUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateProbability",
            (double(DiffusionCaUpdateRule2::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<2> &, double, double, ::CellPtr)) &DiffusionCaUpdateRule2::EvaluateProbability,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX"), py::arg("cell") )
        .def(
            "GetDiffusionParameter",
            (double(DiffusionCaUpdateRule2::*)()) &DiffusionCaUpdateRule2::GetDiffusionParameter,
            " "  )
        .def(
            "SetDiffusionParameter",
            (void(DiffusionCaUpdateRule2::*)(double)) &DiffusionCaUpdateRule2::SetDiffusionParameter,
            " " , py::arg("diffusionParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(DiffusionCaUpdateRule2::*)(::out_stream &)) &DiffusionCaUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
