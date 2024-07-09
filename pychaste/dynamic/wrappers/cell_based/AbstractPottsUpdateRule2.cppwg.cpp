#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPottsUpdateRule.hpp"

#include "AbstractPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPottsUpdateRule<2 > AbstractPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPottsUpdateRule2_Overrides : public AbstractPottsUpdateRule2{
    public:
    using AbstractPottsUpdateRule2::AbstractPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractPottsUpdateRule2,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractPottsUpdateRule2_class(py::module &m){
py::class_<AbstractPottsUpdateRule2 , AbstractPottsUpdateRule2_Overrides , boost::shared_ptr<AbstractPottsUpdateRule2 >  , AbstractUpdateRule<2>  >(m, "AbstractPottsUpdateRule2")
        .def(
            "EvaluateHamiltonianContribution",
            (double(AbstractPottsUpdateRule2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &AbstractPottsUpdateRule2::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractPottsUpdateRule2::*)(::out_stream &)) &AbstractPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
