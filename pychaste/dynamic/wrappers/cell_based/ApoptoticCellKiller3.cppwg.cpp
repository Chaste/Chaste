#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ApoptoticCellKiller.hpp"

#include "ApoptoticCellKiller3.cppwg.hpp"

namespace py = pybind11;
typedef ApoptoticCellKiller<3 > ApoptoticCellKiller3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ApoptoticCellKiller3_Overrides : public ApoptoticCellKiller3{
    public:
    using ApoptoticCellKiller3::ApoptoticCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller3,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ApoptoticCellKiller3,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_ApoptoticCellKiller3_class(py::module &m){
py::class_<ApoptoticCellKiller3 , ApoptoticCellKiller3_Overrides , boost::shared_ptr<ApoptoticCellKiller3 >  , AbstractCellKiller<3>  >(m, "ApoptoticCellKiller3")
        .def(py::init<::AbstractCellPopulation<3> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelSingleCellForApoptosis",
            (void(ApoptoticCellKiller3::*)(::CellPtr)) &ApoptoticCellKiller3::CheckAndLabelSingleCellForApoptosis,
            " " , py::arg("pCell") )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(ApoptoticCellKiller3::*)()) &ApoptoticCellKiller3::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(ApoptoticCellKiller3::*)(::out_stream &)) &ApoptoticCellKiller3::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
