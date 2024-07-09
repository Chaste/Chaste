#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaUpdateRule.hpp"

#include "AbstractCaUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaUpdateRule<3 > AbstractCaUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaUpdateRule3_Overrides : public AbstractCaUpdateRule3{
    public:
    using AbstractCaUpdateRule3::AbstractCaUpdateRule;
    double EvaluateProbability(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::CaBasedCellPopulation<3> & rCellPopulation, double dt, double deltaX, ::CellPtr cell) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaUpdateRule3,
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
            AbstractCaUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaUpdateRule3_class(py::module &m){
py::class_<AbstractCaUpdateRule3 , AbstractCaUpdateRule3_Overrides , boost::shared_ptr<AbstractCaUpdateRule3 >  , AbstractUpdateRule<3>  >(m, "AbstractCaUpdateRule3")
        .def(
            "EvaluateProbability",
            (double(AbstractCaUpdateRule3::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<3> &, double, double, ::CellPtr)) &AbstractCaUpdateRule3::EvaluateProbability,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX"), py::arg("cell") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractCaUpdateRule3::*)(::out_stream &)) &AbstractCaUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
