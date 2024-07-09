#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "IsolatedLabelledCellKiller.hpp"

#include "IsolatedLabelledCellKiller3.cppwg.hpp"

namespace py = pybind11;
typedef IsolatedLabelledCellKiller<3 > IsolatedLabelledCellKiller3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class IsolatedLabelledCellKiller3_Overrides : public IsolatedLabelledCellKiller3{
    public:
    using IsolatedLabelledCellKiller3::IsolatedLabelledCellKiller;
    void CheckAndLabelCellsForApoptosisOrDeath() override {
        PYBIND11_OVERRIDE(
            void,
            IsolatedLabelledCellKiller3,
            CheckAndLabelCellsForApoptosisOrDeath,
            );
    }
    void OutputCellKillerParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            IsolatedLabelledCellKiller3,
            OutputCellKillerParameters,
                    rParamsFile);
    }

};
void register_IsolatedLabelledCellKiller3_class(py::module &m){
py::class_<IsolatedLabelledCellKiller3 , IsolatedLabelledCellKiller3_Overrides , boost::shared_ptr<IsolatedLabelledCellKiller3 >  , AbstractCellKiller<3>  >(m, "IsolatedLabelledCellKiller3")
        .def(py::init<::AbstractCellPopulation<3> * >(), py::arg("pCellPopulation"))
        .def(
            "CheckAndLabelCellsForApoptosisOrDeath",
            (void(IsolatedLabelledCellKiller3::*)()) &IsolatedLabelledCellKiller3::CheckAndLabelCellsForApoptosisOrDeath,
            " "  )
        .def(
            "OutputCellKillerParameters",
            (void(IsolatedLabelledCellKiller3::*)(::out_stream &)) &IsolatedLabelledCellKiller3::OutputCellKillerParameters,
            " " , py::arg("rParamsFile") )
    ;
}
