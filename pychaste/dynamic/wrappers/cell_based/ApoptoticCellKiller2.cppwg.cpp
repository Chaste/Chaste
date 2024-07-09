#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ApoptoticCellKiller.hpp"

#include "ApoptoticCellKiller2.cppwg.hpp"

namespace py = pybind11;
typedef ApoptoticCellKiller<2 > ApoptoticCellKiller2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ApoptoticCellKiller2_Overrides : public ApoptoticCellKiller2{
    public:
    using ApoptoticCellKiller2::ApoptoticCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller2,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_ApoptoticCellKiller2_class(py::module &m){
py::class_<ApoptoticCellKiller2 , ApoptoticCellKiller2_Overrides , boost::shared_ptr<ApoptoticCellKiller2 >  , AbstractCellKiller<2>  >(m, "ApoptoticCellKiller2")
        .def(py::init<::AbstractCellPopulation<2> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelSingleCellForApoptosis",
            (void(ApoptoticCellKiller2::*)(::CellPtr)) &ApoptoticCellKiller2::CheckAndLabelSingleCellForApoptosis,
            " " , py::arg("pCell") )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(ApoptoticCellKiller2::*)()) &ApoptoticCellKiller2::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(ApoptoticCellKiller2::*)(::out_stream &)) &ApoptoticCellKiller2::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
