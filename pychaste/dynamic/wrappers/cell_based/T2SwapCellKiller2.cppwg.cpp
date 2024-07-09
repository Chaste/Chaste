#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "T2SwapCellKiller.hpp"

#include "T2SwapCellKiller2.cppwg.hpp"

namespace py = pybind11;
typedef T2SwapCellKiller<2 > T2SwapCellKiller2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class T2SwapCellKiller2_Overrides : public T2SwapCellKiller2{
    public:
    using T2SwapCellKiller2::T2SwapCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            T2SwapCellKiller2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            T2SwapCellKiller2,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_T2SwapCellKiller2_class(py::module &m){
py::class_<T2SwapCellKiller2 , T2SwapCellKiller2_Overrides , boost::shared_ptr<T2SwapCellKiller2 >  , AbstractCellKiller<2>  >(m, "T2SwapCellKiller2")
        .def(py::init<::AbstractCellPopulation<2> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(T2SwapCellKiller2::*)()) &T2SwapCellKiller2::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(T2SwapCellKiller2::*)(::out_stream &)) &T2SwapCellKiller2::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
