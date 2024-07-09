#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaUpdateRule.hpp"

#include "AbstractCaUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaUpdateRule<2 > AbstractCaUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaUpdateRule2_Overrides : public AbstractCaUpdateRule2{
    public:
    using AbstractCaUpdateRule2::AbstractCaUpdateRule;
    double EvaluateProbability(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::CaBasedCellPopulation<2> & rCellPopulation, double dt, double deltaX, ::CellPtr cell) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaUpdateRule2,
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
            AbstractCaUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaUpdateRule2_class(py::module &m){
py::class_<AbstractCaUpdateRule2 , AbstractCaUpdateRule2_Overrides , boost::shared_ptr<AbstractCaUpdateRule2 >  , AbstractUpdateRule<2>  >(m, "AbstractCaUpdateRule2")
        .def(
            "EvaluateProbability",
            (double(AbstractCaUpdateRule2::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<2> &, double, double, ::CellPtr)) &AbstractCaUpdateRule2::EvaluateProbability,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX"), py::arg("cell") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractCaUpdateRule2::*)(::out_stream &)) &AbstractCaUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
