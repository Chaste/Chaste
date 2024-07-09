#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "T2SwapCellKiller.hpp"

#include "T2SwapCellKiller3.cppwg.hpp"

namespace py = pybind11;
typedef T2SwapCellKiller<3 > T2SwapCellKiller3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class T2SwapCellKiller3_Overrides : public T2SwapCellKiller3{
    public:
    using T2SwapCellKiller3::T2SwapCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            T2SwapCellKiller3,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            T2SwapCellKiller3,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_T2SwapCellKiller3_class(py::module &m){
py::class_<T2SwapCellKiller3 , T2SwapCellKiller3_Overrides , boost::shared_ptr<T2SwapCellKiller3 >  , AbstractCellKiller<3>  >(m, "T2SwapCellKiller3")
        .def(py::init<::AbstractCellPopulation<3> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(T2SwapCellKiller3::*)()) &T2SwapCellKiller3::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(T2SwapCellKiller3::*)(::out_stream &)) &T2SwapCellKiller3::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
