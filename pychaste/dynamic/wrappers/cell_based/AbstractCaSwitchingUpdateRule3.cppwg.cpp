#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaSwitchingUpdateRule.hpp"

#include "AbstractCaSwitchingUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaSwitchingUpdateRule<3 > AbstractCaSwitchingUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCaSwitchingUpdateRule3_Overrides : public AbstractCaSwitchingUpdateRule3{
    public:
    using AbstractCaSwitchingUpdateRule3::AbstractCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<3> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCaSwitchingUpdateRule3,
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
            AbstractCaSwitchingUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaSwitchingUpdateRule3_class(py::module &m){
py::class_<AbstractCaSwitchingUpdateRule3 , AbstractCaSwitchingUpdateRule3_Overrides , boost::shared_ptr<AbstractCaSwitchingUpdateRule3 >  , AbstractUpdateRule<3>  >(m, "AbstractCaSwitchingUpdateRule3")
        .def(
            "EvaluateSwitchingProbability",
            (double(AbstractCaSwitchingUpdateRule3::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<3> &, double, double)) &AbstractCaSwitchingUpdateRule3::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractCaSwitchingUpdateRule3::*)(::out_stream &)) &AbstractCaSwitchingUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
