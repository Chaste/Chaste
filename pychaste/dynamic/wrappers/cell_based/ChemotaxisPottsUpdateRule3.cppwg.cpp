#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"

#include "ChemotaxisPottsUpdateRule3.cppwg.hpp"

namespace py = pybind11;
typedef ChemotaxisPottsUpdateRule<3 > ChemotaxisPottsUpdateRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChemotaxisPottsUpdateRule3_Overrides : public ChemotaxisPottsUpdateRule3{
    public:
    using ChemotaxisPottsUpdateRule3::ChemotaxisPottsUpdateRule;
    double EvaluateHamiltonianContribution(unsigned int currentNodeIndex, unsigned int targetNodeIndex, ::PottsBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ChemotaxisPottsUpdateRule3,
            EvaluateHamiltonianContribution,
                    currentNodeIndex,
        targetNodeIndex,
        rCellPopulation);
    }
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotaxisPottsUpdateRule3,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_ChemotaxisPottsUpdateRule3_class(py::module &m){
py::class_<ChemotaxisPottsUpdateRule3 , ChemotaxisPottsUpdateRule3_Overrides , boost::shared_ptr<ChemotaxisPottsUpdateRule3 >  , AbstractPottsUpdateRule<3>  >(m, "ChemotaxisPottsUpdateRule3")
        .def(py::init< >())
        .def(
            "EvaluateHamiltonianContribution",
            (double(ChemotaxisPottsUpdateRule3::*)(unsigned int, unsigned int, ::PottsBasedCellPopulation<3> &)) &ChemotaxisPottsUpdateRule3::EvaluateHamiltonianContribution,
            " " , py::arg("currentNodeIndex"), py::arg("targetNodeIndex"), py::arg("rCellPopulation") )
        .def(
            "OutputUpdateRuleParameters",
            (void(ChemotaxisPottsUpdateRule3::*)(::out_stream &)) &ChemotaxisPottsUpdateRule3::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
