#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "OffLatticeSimulation.hpp"

#include "OffLatticeSimulation2_2.cppwg.hpp"

namespace py = pybind11;
typedef OffLatticeSimulation<2,2 > OffLatticeSimulation2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class OffLatticeSimulation2_2_Overrides : public OffLatticeSimulation2_2{
    public:
    using OffLatticeSimulation2_2::OffLatticeSimulation;
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation2_2,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation2_2,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation2_2,
            UpdateCellLocationsAndTopology,
            );
    }
    void SetupSolve() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation2_2,
            SetupSolve,
            );
    }
    void WriteVisualizerSetupFile() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation2_2,
            WriteVisualizerSetupFile,
            );
    }

};
void register_OffLatticeSimulation2_2_class(py::module &m){
py::class_<OffLatticeSimulation2_2 , OffLatticeSimulation2_2_Overrides , boost::shared_ptr<OffLatticeSimulation2_2 >  , AbstractCellBasedSimulation<2, 2>  >(m, "OffLatticeSimulation2_2")
        .def(py::init<::AbstractCellPopulation<2> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "AddForce",
            (void(OffLatticeSimulation2_2::*)(::boost::shared_ptr<AbstractForce<2, 2>>)) &OffLatticeSimulation2_2::AddForce,
            " " , py::arg("pForce") )
        .def(
            "RemoveAllForces",
            (void(OffLatticeSimulation2_2::*)()) &OffLatticeSimulation2_2::RemoveAllForces,
            " "  )
        .def(
            "AddCellPopulationBoundaryCondition",
            (void(OffLatticeSimulation2_2::*)(::boost::shared_ptr<AbstractCellPopulationBoundaryCondition<2, 2>>)) &OffLatticeSimulation2_2::AddCellPopulationBoundaryCondition,
            " " , py::arg("pBoundaryCondition") )
        .def(
            "RemoveAllCellPopulationBoundaryConditions",
            (void(OffLatticeSimulation2_2::*)()) &OffLatticeSimulation2_2::RemoveAllCellPopulationBoundaryConditions,
            " "  )
        .def(
            "SetNumericalMethod",
            (void(OffLatticeSimulation2_2::*)(::boost::shared_ptr<AbstractNumericalMethod<2, 2>>)) &OffLatticeSimulation2_2::SetNumericalMethod,
            " " , py::arg("pNumericalMethod") )
        .def(
            "GetNumericalMethod",
            (::boost::shared_ptr<AbstractNumericalMethod<2, 2>> const(OffLatticeSimulation2_2::*)() const ) &OffLatticeSimulation2_2::GetNumericalMethod,
            " "  )
        .def(
            "OutputAdditionalSimulationSetup",
            (void(OffLatticeSimulation2_2::*)(::out_stream &)) &OffLatticeSimulation2_2::OutputAdditionalSimulationSetup,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationParameters",
            (void(OffLatticeSimulation2_2::*)(::out_stream &)) &OffLatticeSimulation2_2::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "rGetForceCollection",
            (::std::vector<boost::shared_ptr<AbstractForce<2, 2>>> const &(OffLatticeSimulation2_2::*)() const ) &OffLatticeSimulation2_2::rGetForceCollection,
            " "  , py::return_value_policy::reference_internal)
    ;
}
