#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedSimulation.hpp"

#include "AbstractCellBasedSimulation3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedSimulation<3,3 > AbstractCellBasedSimulation3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractCellBasedSimulation3_3_Overrides : public AbstractCellBasedSimulation3_3{
    public:
    using AbstractCellBasedSimulation3_3::AbstractCellBasedSimulation;
    void OutputSimulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation3_3,
            OutputSimulationParameters,
                    rParamsFile);
    }
    void WriteVisualizerSetupFile() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation3_3,
            WriteVisualizerSetupFile,
            );
    }
    unsigned int DoCellBirth() override {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractCellBasedSimulation3_3,
            DoCellBirth,
            );
    }
    void SetupSolve() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation3_3,
            SetupSolve,
            );
    }
    bool StoppingEventHasOccurred() override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellBasedSimulation3_3,
            StoppingEventHasOccurred,
            );
    }
    void UpdateCellPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedSimulation3_3,
            UpdateCellPopulation,
            );
    }
    void UpdateCellLocationsAndTopology() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation3_3,
            UpdateCellLocationsAndTopology,
            );
    }
    void OutputAdditionalSimulationSetup(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellBasedSimulation3_3,
            OutputAdditionalSimulationSetup,
                    rParamsFile);
    }

};
void register_AbstractCellBasedSimulation3_3_class(py::module &m){
py::class_<AbstractCellBasedSimulation3_3 , AbstractCellBasedSimulation3_3_Overrides , boost::shared_ptr<AbstractCellBasedSimulation3_3 >   >(m, "AbstractCellBasedSimulation3_3")
        .def(py::init<::AbstractCellPopulation<3> &, bool, bool >(), py::arg("rCellPopulation"), py::arg("deleteCellPopulationInDestructor") = false, py::arg("initialiseCells") = true)
        .def(
            "GetNodeLocation",
            (::std::vector<double>(AbstractCellBasedSimulation3_3::*)(unsigned int const &)) &AbstractCellBasedSimulation3_3::GetNodeLocation,
            " " , py::arg("rNodeIndex") )
        .def(
            "GetDt",
            (double(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetDt,
            " "  )
        .def(
            "GetNumBirths",
            (unsigned int(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetNumBirths,
            " "  )
        .def(
            "GetNumDeaths",
            (unsigned int(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetNumDeaths,
            " "  )
        .def(
            "GetOutputDirectory",
            (::std::string(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetOutputDirectory,
            " "  )
        .def(
            "SetDt",
            (void(AbstractCellBasedSimulation3_3::*)(double)) &AbstractCellBasedSimulation3_3::SetDt,
            " " , py::arg("dt") )
        .def(
            "SetEndTime",
            (void(AbstractCellBasedSimulation3_3::*)(double)) &AbstractCellBasedSimulation3_3::SetEndTime,
            " " , py::arg("endTime") )
        .def(
            "SetOutputDirectory",
            (void(AbstractCellBasedSimulation3_3::*)(::std::string)) &AbstractCellBasedSimulation3_3::SetOutputDirectory,
            " " , py::arg("outputDirectory") )
        .def(
            "SetSamplingTimestepMultiple",
            (void(AbstractCellBasedSimulation3_3::*)(unsigned int)) &AbstractCellBasedSimulation3_3::SetSamplingTimestepMultiple,
            " " , py::arg("samplingTimestepMultiple") )
        .def(
            "SetUpdatingTimestepMultiple",
            (void(AbstractCellBasedSimulation3_3::*)(unsigned int)) &AbstractCellBasedSimulation3_3::SetUpdatingTimestepMultiple,
            " " , py::arg("updatingTimestepMultiple") )
        .def(
            "SetNoBirth",
            (void(AbstractCellBasedSimulation3_3::*)(bool)) &AbstractCellBasedSimulation3_3::SetNoBirth,
            " " , py::arg("noBirth") )
        .def(
            "SetUpdateCellPopulationRule",
            (void(AbstractCellBasedSimulation3_3::*)(bool)) &AbstractCellBasedSimulation3_3::SetUpdateCellPopulationRule,
            " " , py::arg("updateCellPopulation") )
        .def(
            "GetUpdateCellPopulationRule",
            (bool(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetUpdateCellPopulationRule,
            " "  )
        .def(
            "AddCellKiller",
            (void(AbstractCellBasedSimulation3_3::*)(::boost::shared_ptr<AbstractCellKiller<3>>)) &AbstractCellBasedSimulation3_3::AddCellKiller,
            " " , py::arg("pCellKiller") )
        .def(
            "RemoveAllCellKillers",
            (void(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::RemoveAllCellKillers,
            " "  )
        .def(
            "AddSimulationModifier",
            (void(AbstractCellBasedSimulation3_3::*)(::boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>)) &AbstractCellBasedSimulation3_3::AddSimulationModifier,
            " " , py::arg("pSimulationModifier") )
        .def(
            "AddTopologyUpdateSimulationModifier",
            (void(AbstractCellBasedSimulation3_3::*)(::boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>)) &AbstractCellBasedSimulation3_3::AddTopologyUpdateSimulationModifier,
            " " , py::arg("pSimulationModifier") )
        .def(
            "GetTopologyUpdateSimulationModifiers",
            (::std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3>>> *(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetTopologyUpdateSimulationModifiers,
            " "  , py::return_value_policy::reference)
        .def(
            "Solve",
            (void(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::Solve,
            " "  )
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<3> &(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetCellPopulation",
            (::AbstractCellPopulation<3> const &(AbstractCellBasedSimulation3_3::*)() const ) &AbstractCellBasedSimulation3_3::rGetCellPopulation,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetOutputDivisionLocations",
            (bool(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetOutputDivisionLocations,
            " "  )
        .def(
            "SetOutputDivisionLocations",
            (void(AbstractCellBasedSimulation3_3::*)(bool)) &AbstractCellBasedSimulation3_3::SetOutputDivisionLocations,
            " " , py::arg("outputDivisionLocations") )
        .def(
            "GetOutputCellVelocities",
            (bool(AbstractCellBasedSimulation3_3::*)()) &AbstractCellBasedSimulation3_3::GetOutputCellVelocities,
            " "  )
        .def(
            "SetOutputCellVelocities",
            (void(AbstractCellBasedSimulation3_3::*)(bool)) &AbstractCellBasedSimulation3_3::SetOutputCellVelocities,
            " " , py::arg("outputCellVelocities") )
        .def(
            "OutputSimulationParameters",
            (void(AbstractCellBasedSimulation3_3::*)(::out_stream &)) &AbstractCellBasedSimulation3_3::OutputSimulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
