#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOdeBasedPhaseBasedCellCycleModel.hpp"

#include "AbstractOdeBasedPhaseBasedCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeBasedPhaseBasedCellCycleModel AbstractOdeBasedPhaseBasedCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeBasedPhaseBasedCellCycleModel_Overrides : public AbstractOdeBasedPhaseBasedCellCycleModel{
    public:
    using AbstractOdeBasedPhaseBasedCellCycleModel::AbstractOdeBasedPhaseBasedCellCycleModel;
    void UpdateCellCyclePhase() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedPhaseBasedCellCycleModel,
            UpdateCellCyclePhase,
            );
    }
    void SetBirthTime(double birthTime) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedPhaseBasedCellCycleModel,
            SetBirthTime,
                    birthTime);
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedPhaseBasedCellCycleModel,
            ResetForDivision,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedPhaseBasedCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AbstractOdeBasedPhaseBasedCellCycleModel_class(py::module &m){
py::class_<AbstractOdeBasedPhaseBasedCellCycleModel , AbstractOdeBasedPhaseBasedCellCycleModel_Overrides , boost::shared_ptr<AbstractOdeBasedPhaseBasedCellCycleModel >  , AbstractPhaseBasedCellCycleModel  >(m, "AbstractOdeBasedPhaseBasedCellCycleModel")
        .def(
            "UpdateCellCyclePhase",
            (void(AbstractOdeBasedPhaseBasedCellCycleModel::*)()) &AbstractOdeBasedPhaseBasedCellCycleModel::UpdateCellCyclePhase,
            " "  )
        .def(
            "GetOdeStopTime",
            (double(AbstractOdeBasedPhaseBasedCellCycleModel::*)()) &AbstractOdeBasedPhaseBasedCellCycleModel::GetOdeStopTime,
            " "  )
        .def(
            "SetBirthTime",
            (void(AbstractOdeBasedPhaseBasedCellCycleModel::*)(double)) &AbstractOdeBasedPhaseBasedCellCycleModel::SetBirthTime,
            " " , py::arg("birthTime") )
        .def(
            "ResetForDivision",
            (void(AbstractOdeBasedPhaseBasedCellCycleModel::*)()) &AbstractOdeBasedPhaseBasedCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(AbstractOdeBasedPhaseBasedCellCycleModel::*)(::out_stream &)) &AbstractOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
