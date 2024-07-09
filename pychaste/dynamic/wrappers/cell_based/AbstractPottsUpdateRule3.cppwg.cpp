#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPottsUpdateRule.hpp"

#include "AbstractPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPottsUpdateRule<3 > AbstractPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPottsUpdateRule3_Overrides : public AbstractPottsUpdateRule3{
    public:
    using AbstractPottsUpdateRule3::AbstractPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractPottsUpdateRule3_class(py::module &m){
py::class_<AbstractPottsUpdateRule3 , AbstractPottsUpdateRule3_Overrides , boost::shared_ptr<AbstractPottsUpdateRule3 >  , AbstractUpdateRule<3>  >(m, "AbstractPottsUpdateRule3")
        .def(
            "EvaluateHamiltonianContribution",
            (double(AbstractPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &AbstractPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractPottsUpdateRule3::*)(::out_stream &)) &AbstractPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
