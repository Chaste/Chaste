#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"

#include "RandomCaSwitchingUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef RandomCaSwitchingUpdateRule<3 > RandomCaSwitchingUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RandomCaSwitchingUpdateRule3_Overrides : public RandomCaSwitchingUpdateRule3{
    public:
    using RandomCaSwitchingUpdateRule3::RandomCaSwitchingUpdateRule;
    double EvaluateSwitchingProbability(unsigned int currentNodeIndex, unsigned int neighbourNodeIndex, ::CaBasedCellPopulation<3> & rCellPopulation, double dt, double deltaX) override {
        PYBIND11_OVERRIDE(
            double,
            RandomCaSwitchingUpdateRule3,
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
            RandomCaSwitchingUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_RandomCaSwitchingUpdateRule3_class(py::module &m){
py::class_<RandomCaSwitchingUpdateRule3 , RandomCaSwitchingUpdateRule3_Overrides , boost::shared_ptr<RandomCaSwitchingUpdateRule3 >  , AbstractCaSwitchingUpdateRule<3>  >(m, "RandomCaSwitchingUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateSwitchingProbability",
            (double(RandomCaSwitchingUpdateRule3::*)(unsigned int, unsigned int, ::CaBasedCellPopulation<3> &, double, double)) &RandomCaSwitchingUpdateRule3::EvaluateSwitchingProbability,
            " " , py::arg("currentNodeIndex"), py::arg("neighbourNodeIndex"), py::arg("rCellPopulation"), py::arg("dt"), py::arg("deltaX") )
        .def(
            "GetSwitchingParameter",
            (double(RandomCaSwitchingUpdateRule3::*)()) &RandomCaSwitchingUpdateRule3::GetSwitchingParameter,
            " "  )
        .def(
            "SetSwitchingParameter",
            (void(RandomCaSwitchingUpdateRule3::*)(double)) &RandomCaSwitchingUpdateRule3::SetSwitchingParameter,
            " " , py::arg("switchingParameter") )
        .def(
            "OutputUpdateRuleParameters",
            (void(RandomCaSwitchingUpdateRule3::*)(::out_stream &)) &RandomCaSwitchingUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
