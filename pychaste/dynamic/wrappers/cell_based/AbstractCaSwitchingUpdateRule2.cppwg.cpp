#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaSwitchingUpdateRule.hpp"

#include "AbstractCaSwitchingUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaSwitchingUpdateRule<2 > AbstractCaSwitchingUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaSwitchingUpdateRule2_Overrides : public AbstractCaSwitchingUpdateRule2{
    public:
    using AbstractCaSwitchingUpdateRule2::AbstractCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<2> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaSwitchingUpdateRule2,
            EvaluateSwitchingProbability,
                    currentNodeIndex,
        neighbourNodeIndex,
        rCellPopulation,
        dt,
        deltaX);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCaSwitchingUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaSwitchingUpdateRule2_class(py::module &m){
py::class_<AbstractCaSwitchingUpdateRule2 , AbstractCaSwitchingUpdateRule2_Overrides , boost::shared_ptr<AbstractCaSwitchingUpdateRule2 >  , AbstractUpdateRule<2>  >(m, "AbstractCaSwitchingUpdateRule2")
        .def(
            "EvaluateSwitchingProbability",
            (double(AbstractCaSwitchingUpdateRule2::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<2> &, double, double)) &AbstractCaSwitchingUpdateRule2::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractCaSwitchingUpdateRule2::*)(::out_stream &)) &AbstractCaSwitchingUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
