#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RandomCellKiller.hpp"

#include "RandomCellKiller2.cppwg.hpp"

namespace py = pybind11;
typedef RandomCellKiller<2 > RandomCellKiller2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RandomCellKiller2_Overrides : public RandomCellKiller2{
    public:
    using RandomCellKiller2::RandomCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            RandomCellKiller2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            RandomCellKiller2,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_RandomCellKiller2_class(py::module &m){
py::class_<RandomCellKiller2 , RandomCellKiller2_Overrides , boost::shared_ptr<RandomCellKiller2 >  , AbstractCellKiller<2>  >(m, "RandomCellKiller2")
        .def(py::init<::AbstractCellPopulation<2> *, double >(), py::arg("pCellPopulation"), py::arg("probabilityOfDeathInAnHour"))
        .def(
            "GetDeathProbabilityInAnHour",
            (double(RandomCellKiller2::*)() const ) &RandomCellKiller2::GetDeathProbabilityInAnHour,
            " "  )
        .def(
            "CheckAndLabelSingleCellForApoptosis",
            (void(RandomCellKiller2::*)(::CellPtr)) &RandomCellKiller2::CheckAndLabelSingleCellForApoptosis,
            " " , py::arg("pCell") )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(RandomCellKiller2::*)()) &RandomCellKiller2::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(RandomCellKiller2::*)(::out_stream &)) &RandomCellKiller2::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
