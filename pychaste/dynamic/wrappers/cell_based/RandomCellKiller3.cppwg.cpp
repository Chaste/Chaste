#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomCellKiller.hpp"

#include "RandomCellKiller3.cppwg.hpp"

namespace py = pybind11;
typedef RandomCellKiller<3 > RandomCellKiller3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RandomCellKiller3_Overrides : public RandomCellKiller3{
    public:
    using RandomCellKiller3::RandomCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            RandomCellKiller3,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            RandomCellKiller3,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_RandomCellKiller3_class(py::module &m){
py::class_<RandomCellKiller3 , RandomCellKiller3_Overrides , boost::shared_ptr<RandomCellKiller3 >  , AbstractCellKiller<3>  >(m, "RandomCellKiller3")
        .def(py::init<::AbstractCellPopulation<3> *, double >(), py::arg("pCellPopulation"), py::arg("probabilityOfDeathInAnHour"))
        .def(
            "GetDeathProbabilityInAnHour",
            (double(RandomCellKiller3::*)() const ) &RandomCellKiller3::GetDeathProbabilityInAnHour,
            " "  )
        .def(
            "CheckAndLabelSingleCellForApoptosis",
            (void(RandomCellKiller3::*)(::CellPtr)) &RandomCellKiller3::CheckAndLabelSingleCellForApoptosis,
            " " , py::arg("pCell") )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(RandomCellKiller3::*)()) &RandomCellKiller3::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(RandomCellKiller3::*)(::out_stream &)) &RandomCellKiller3::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
