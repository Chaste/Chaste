#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"

#include "RandomCaSwitchingUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef RandomCaSwitchingUpdateRule<2 > RandomCaSwitchingUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RandomCaSwitchingUpdateRule2_Overrides : public RandomCaSwitchingUpdateRule2{
    public:
    using RandomCaSwitchingUpdateRule2::RandomCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<2> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE(
            double,
            RandomCaSwitchingUpdateRule2,
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
            RandomCaSwitchingUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_RandomCaSwitchingUpdateRule2_class(py::module &m){
py::class_<RandomCaSwitchingUpdateRule2 , RandomCaSwitchingUpdateRule2_Overrides , boost::shared_ptr<RandomCaSwitchingUpdateRule2 >  , AbstractCaSwitchingUpdateRule<2>  >(m, "RandomCaSwitchingUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateSwitchingProbability",
            (double(RandomCaSwitchingUpdateRule2::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<2> &, double, double)) &RandomCaSwitchingUpdateRule2::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "GetSwitchingParameter",
            (double(RandomCaSwitchingUpdateRule2::*)()) &RandomCaSwitchingUpdateRule2::GetSwitchingParameter,
            " "  )
        .def(
            "SetSwitchingParameter",
            (void(RandomCaSwitchingUpdateRule2::*)(double)) &RandomCaSwitchingUpdateRule2::SetSwitchingParameter,
            " " , py::arg("switchingParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(RandomCaSwitchingUpdateRule2::*)(::out_stream &)) &RandomCaSwitchingUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
