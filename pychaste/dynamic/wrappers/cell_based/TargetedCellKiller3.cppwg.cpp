#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TargetedCellKiller.hpp"

#include "TargetedCellKiller3.cppwg.hpp"

namespace py = pybind11;
typedef TargetedCellKiller<3 > TargetedCellKiller3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class TargetedCellKiller3_Overrides : public TargetedCellKiller3{
    public:
    using TargetedCellKiller3::TargetedCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            TargetedCellKiller3,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            TargetedCellKiller3,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_TargetedCellKiller3_class(py::module &m){
py::class_<TargetedCellKiller3 , TargetedCellKiller3_Overrides , boost::shared_ptr<TargetedCellKiller3 >  , AbstractCellKiller<3>  >(m, "TargetedCellKiller3")
        .def(py::init<::AbstractCellPopulation<3> *, unsigned int, bool >(), py::arg("pCellPopulation"), py::arg("targetedIndex"), py::arg("bloodLust") = true)
        .def(
            "GetTargetIndex",
            (unsigned int(TargetedCellKiller3::*)() const ) &TargetedCellKiller3::GetTargetIndex,
            " "  )
        .def(
            "GetBloodLust",
            (unsigned int(TargetedCellKiller3::*)() const ) &TargetedCellKiller3::GetBloodLust,
            " "  )
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(TargetedCellKiller3::*)()) &TargetedCellKiller3::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(TargetedCellKiller3::*)(::out_stream &)) &TargetedCellKiller3::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
