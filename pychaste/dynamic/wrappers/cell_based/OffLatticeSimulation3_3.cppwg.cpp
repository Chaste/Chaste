#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "OffLatticeSimulation.hpp"

#include "OffLatticeSimulation3_3.cppwg.hpp"

namespace py = pybind11;
typedef OffLatticeSimulation<3,3 > OffLatticeSimulation3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class OffLatticeSimulation3_3_Overrides : public OffLatticeSimulation3_3{
    public:
    using OffLatticeSimulation3_3::OffLatticeSimulation;
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation3_3,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation3_3,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation3_3,
            UpdateCellLocationsAndTopology,
            );
    }
    void SetupSolve() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation3_3,
            SetupSolve,
            );
    }
    void WriteVisualizerSetupFile() override {
        PYBIND11_OVERRIDE(
            void,
            OffLatticeSimulation3_3,
            WriteVisualizerSetupFile,
            );
    }

};
void register_OffLatticeSimulation3_3_class(py::module &m){
py::class_<OffLatticeSimulation3_3 , OffLatticeSimulation3_3_Overrides , boost::shared_ptr<OffLatticeSimulation3_3 >  , AbstractCellBasedSimulation<3, 3>  >(m, "OffLatticeSimulation3_3")
        .def(py::init<::AbstractCellPopulation<3> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "AddForce",
            (void(OffLatticeSimulation3_3::*)(::boost::shared_ptr<AbstractForce<3, 3>>)) &OffLatticeSimulation3_3::AddForce,
            " " , py::arg("pForce") )
        .def(
            "RemoveAllForces",
            (void(OffLatticeSimulation3_3::*)()) &OffLatticeSimulation3_3::RemoveAllForces,
            " "  )
        .def(
            "AddCellPopulationBoundaryCondition",
            (void(OffLatticeSimulation3_3::*)(::boost::shared_ptr<AbstractCellPopulationBoundaryCondition<3, 3>>)) &OffLatticeSimulation3_3::AddCellPopulationBoundaryCondition,
            " " , py::arg("pBoundaryCondition") )
        .def(
            "RemoveAllCellPopulationBoundaryConditions",
            (void(OffLatticeSimulation3_3::*)()) &OffLatticeSimulation3_3::RemoveAllCellPopulationBoundaryConditions,
            " "  )
        .def(
            "SetNumericalMethod",
            (void(OffLatticeSimulation3_3::*)(::boost::shared_ptr<AbstractNumericalMethod<3, 3>>)) &OffLatticeSimulation3_3::SetNumericalMethod,
            " " , py::arg("pNumericalMethod") )
        .def(
            "GetNumericalMethod",
            (::boost::shared_ptr<AbstractNumericalMethod<3, 3>> const(OffLatticeSimulation3_3::*)() const ) &OffLatticeSimulation3_3::GetNumericalMethod,
            " "  )
        .def(
            "OutputAdditionalSimulationSetup",
            (void(OffLatticeSimulation3_3::*)(::out_stream &)) &OffLatticeSimulation3_3::OutputAdditionalSimulationSetup,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputSimulationParameters",
            (void(OffLatticeSimulation3_3::*)(::out_stream &)) &OffLatticeSimulation3_3::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "rGetForceCollection",
            (::std::vector<boost::shared_ptr<AbstractForce<3, 3>>> const &(OffLatticeSimulation3_3::*)() const ) &OffLatticeSimulation3_3::rGetForceCollection,
            " "  , py::return_value_policy::reference_internal)
    ;
}
