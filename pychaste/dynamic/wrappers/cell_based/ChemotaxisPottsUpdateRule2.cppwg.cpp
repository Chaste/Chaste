#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"

#include "ChemotaxisPottsUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef ChemotaxisPottsUpdateRule<2 > ChemotaxisPottsUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChemotaxisPottsUpdateRule2_Overrides : public ChemotaxisPottsUpdateRule2{
    public:
    using ChemotaxisPottsUpdateRule2::ChemotaxisPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ChemotaxisPottsUpdateRule2,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotaxisPottsUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_ChemotaxisPottsUpdateRule2_class(py::module &m){
py::class_<ChemotaxisPottsUpdateRule2 , ChemotaxisPottsUpdateRule2_Overrides , boost::shared_ptr<ChemotaxisPottsUpdateRule2 >  , AbstractPottsUpdateRule<2>  >(m, "ChemotaxisPottsUpdateRule2")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(ChemotaxisPottsUpdateRule2::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<2> &)) &ChemotaxisPottsUpdateRule2::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "OutputUpdateRuleParameters",
            (void(ChemotaxisPottsUpdateRule2::*)(::out_stream &)) &ChemotaxisPottsUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
