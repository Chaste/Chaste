#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TargetedCellKiller.hpp"

#include "TargetedCellKiller2.cppwg.hpp"

namespace py = pybind11;
typedef TargetedCellKiller<2 > TargetedCellKiller2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TargetedCellKiller2_Overrides : public TargetedCellKiller2{
    public:
    using TargetedCellKiller2::TargetedCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            TargetedCellKiller2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            TargetedCellKiller2,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_TargetedCellKiller2_class(py::module &m){
py::class_<TargetedCellKiller2 , TargetedCellKiller2_Overrides , boost::shared_ptr<TargetedCellKiller2 >  , AbstractCellKiller<2>  >(m, "TargetedCellKiller2")
        .def(py::init<::AbstractCellPopulation<2> *, unsigned int, bool >(), py::arg("pCellPopulation"), py::arg("targetedIndex"), py::arg("bloodLust") = true)
        .def(
            "GetTargetIndex",
            (unsigned int(TargetedCellKiller2::*)() const ) &TargetedCellKiller2::GetTargetIndex,
            " "  )
        .def(
            "GetBloodLust",
            (unsigned int(TargetedCellKiller2::*)() const ) &TargetedCellKiller2::GetBloodLust,
            " "  )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(TargetedCellKiller2::*)()) &TargetedCellKiller2::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(TargetedCellKiller2::*)(::out_stream &)) &TargetedCellKiller2::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
