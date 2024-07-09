#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "OnLatticeSimulation.hpp"

#include "OnLatticeSimulation2.cppwg.hpp"

namespace py = pybind11;
typedef OnLatticeSimulation<2 > OnLatticeSimulation2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class OnLatticeSimulation2_Overrides : public OnLatticeSimulation2{
    public:
    using OnLatticeSimulation2::OnLatticeSimulation;
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void UpdateCellPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            UpdateCellPopulation,
            );
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation2,
            UpdateCellLocationsAndTopology,
            );
    }

};
void register_OnLatticeSimulation2_class(py::module &m){
py::class_<OnLatticeSimulation2 , OnLatticeSimulation2_Overrides , boost::shared_ptr<OnLatticeSimulation2 >  , AbstractCellBasedSimulation<2, 2>  >(m, "OnLatticeSimulation2")
        .def(py::init<::AbstractCellPopulation<2> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "AddUpdateRule",
            (void(OnLatticeSimulation2::*)(::boost::shared_ptr<AbstractUpdateRule<2>>)) &OnLatticeSimulation2::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "RemoveAllUpdateRules",
            (void(OnLatticeSimulation2::*)()) &OnLatticeSimulation2::RemoveAllUpdateRules,
            " "  )
        .def(
            "OutputAdditionalSimulationSetup",
            (void(OnLatticeSimulation2::*)(::out_stream &)) &OnLatticeSimulation2::OutputAdditionalSimulationSetup,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationParameters",
            (void(OnLatticeSimulation2::*)(::out_stream &)) &OnLatticeSimulation2::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
