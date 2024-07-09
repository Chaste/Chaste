#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DivisionBiasTrackingModifier.hpp"

#include "DivisionBiasTrackingModifier2.cppwg.hpp"

namespace py = pybind11;
typedef DivisionBiasTrackingModifier<2 > DivisionBiasTrackingModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DivisionBiasTrackingModifier2_Overrides : public DivisionBiasTrackingModifier2{
    public:
    using DivisionBiasTrackingModifier2::DivisionBiasTrackingModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DivisionBiasTrackingModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_DivisionBiasTrackingModifier2_class(py::module &m){
py::class_<DivisionBiasTrackingModifier2 , DivisionBiasTrackingModifier2_Overrides , boost::shared_ptr<DivisionBiasTrackingModifier2 >  , AbstractCellBasedSimulationModifier<2>  >(m, "DivisionBiasTrackingModifier2")
        .def(py::init<::boost::numeric::ublas::c_vector<double, 2> >(), py::arg("divisionBiasVector"))
        .def(
            "rGetDivisionBiasVector",
            (::boost::numeric::ublas::c_vector<double, 2> const &(DivisionBiasTrackingModifier2::*)() const ) &DivisionBiasTrackingModifier2::rGetDivisionBiasVector,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(DivisionBiasTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DivisionBiasTrackingModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(DivisionBiasTrackingModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &DivisionBiasTrackingModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateCellData",
            (void(DivisionBiasTrackingModifier2::*)(::AbstractCellPopulation<2> &)) &DivisionBiasTrackingModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(DivisionBiasTrackingModifier2::*)(::out_stream &)) &DivisionBiasTrackingModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
