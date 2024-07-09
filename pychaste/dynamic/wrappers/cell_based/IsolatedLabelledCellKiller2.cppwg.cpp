#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "IsolatedLabelledCellKiller.hpp"

#include "IsolatedLabelledCellKiller2.cppwg.hpp"

namespace py = pybind11;
typedef IsolatedLabelledCellKiller<2 > IsolatedLabelledCellKiller2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class IsolatedLabelledCellKiller2_Overrides : public IsolatedLabelledCellKiller2{
    public:
    using IsolatedLabelledCellKiller2::IsolatedLabelledCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            IsolatedLabelledCellKiller2,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            IsolatedLabelledCellKiller2,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_IsolatedLabelledCellKiller2_class(py::module &m){
py::class_<IsolatedLabelledCellKiller2 , IsolatedLabelledCellKiller2_Overrides , boost::shared_ptr<IsolatedLabelledCellKiller2 >  , AbstractCellKiller<2>  >(m, "IsolatedLabelledCellKiller2")
        .def(py::init<::AbstractCellPopulation<2> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(IsolatedLabelledCellKiller2::*)()) &IsolatedLabelledCellKiller2::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(IsolatedLabelledCellKiller2::*)(::out_stream &)) &IsolatedLabelledCellKiller2::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
