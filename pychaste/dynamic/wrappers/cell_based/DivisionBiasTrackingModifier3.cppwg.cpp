#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DivisionBiasTrackingModifier.hpp"

#include "DivisionBiasTrackingModifier3.cppwg.hpp"

namespace py = pybind11;
typedef DivisionBiasTrackingModifier<3 > DivisionBiasTrackingModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DivisionBiasTrackingModifier3_Overrides : public DivisionBiasTrackingModifier3{
    public:
    using DivisionBiasTrackingModifier3::DivisionBiasTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DivisionBiasTrackingModifier3_class(py::module &m){
py::class_<DivisionBiasTrackingModifier3 , DivisionBiasTrackingModifier3_Overrides , boost::shared_ptr<DivisionBiasTrackingModifier3 >  , AbstractCellBasedSimulationModifier<3>  >(m, "DivisionBiasTrackingModifier3")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 3> >(), py::arg("divisionBiasVector"))
        .def(
            "rGetDivisionBiasVector",
            (::boost::numeric::ublas::c_vector<double, 3> const &(DivisionBiasTrackingModifier3::*)() const ) &DivisionBiasTrackingModifier3::rGetDivisionBiasVector,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DivisionBiasTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DivisionBiasTrackingModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DivisionBiasTrackingModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &DivisionBiasTrackingModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DivisionBiasTrackingModifier3::*)(::AbstractCellPopulation<3> &)) &DivisionBiasTrackingModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DivisionBiasTrackingModifier3::*)(::out_stream &)) &DivisionBiasTrackingModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
