#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "OnLatticeSimulation.hpp"

#include "OnLatticeSimulation3.cppwg.hpp"

namespace py = pybind11;
typedef OnLatticeSimulation<3 > OnLatticeSimulation3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class OnLatticeSimulation3_Overrides : public OnLatticeSimulation3{
    public:
    using OnLatticeSimulation3::OnLatticeSimulation;
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation3,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation3,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void UpdateCellPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation3,
            UpdateCellPopulation,
            );
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE(
            void,
            OnLatticeSimulation3,
            UpdateCellLocationsAndTopology,
            );
    }

};
void register_OnLatticeSimulation3_class(py::module &m){
py::class_<OnLatticeSimulation3 , OnLatticeSimulation3_Overrides , boost::shared_ptr<OnLatticeSimulation3 >  , AbstractCellBasedSimulation<3, 3>  >(m, "OnLatticeSimulation3")
        .def(py::init<::AbstractCellPopulation<3> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "AddUpdateRule",
            (void(OnLatticeSimulation3::*)(::boost::shared_ptr<AbstractUpdateRule<3>>)) &OnLatticeSimulation3::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "RemoveAllUpdateRules",
            (void(OnLatticeSimulation3::*)()) &OnLatticeSimulation3::RemoveAllUpdateRules,
            " "  )
        .def(
            "OutputAdditionalSimulationSetup",
            (void(OnLatticeSimulation3::*)(::out_stream &)) &OnLatticeSimulation3::OutputAdditionalSimulationSetup,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationParameters",
            (void(OnLatticeSimulation3::*)(::out_stream &)) &OnLatticeSimulation3::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
