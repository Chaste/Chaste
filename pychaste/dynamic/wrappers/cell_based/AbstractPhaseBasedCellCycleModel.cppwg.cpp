#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"

#include "AbstractPhaseBasedCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPhaseBasedCellCycleModel AbstractPhaseBasedCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPhaseBasedCellCycleModel_Overrides : public AbstractPhaseBasedCellCycleModel{
    public:
    using AbstractPhaseBasedCellCycleModel::AbstractPhaseBasedCellCycleModel;
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractPhaseBasedCellCycleModel,
            ReadyToDivide,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPhaseBasedCellCycleModel,
            ResetForDivision,
            );
    }
    void UpdateCellCyclePhase() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPhaseBasedCellCycleModel,
            UpdateCellCyclePhase,
            );
    }
    double GetG1Duration() const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetG1Duration,
            );
    }
    double GetSDuration() const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetSDuration,
            );
    }
    double GetG2Duration() const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetG2Duration,
            );
    }
    double GetMDuration() const  override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetMDuration,
            );
    }
    void SetStemCellG1Duration(double stemCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPhaseBasedCellCycleModel,
            SetStemCellG1Duration,
                    stemCellG1Duration);
    }
    void SetTransitCellG1Duration(double transitCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractPhaseBasedCellCycleModel,
            SetTransitCellG1Duration,
                    transitCellG1Duration);
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractPhaseBasedCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPhaseBasedCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AbstractPhaseBasedCellCycleModel_class(py::module &m){
py::class_<AbstractPhaseBasedCellCycleModel , AbstractPhaseBasedCellCycleModel_Overrides , boost::shared_ptr<AbstractPhaseBasedCellCycleModel >  , AbstractCellCycleModel  >(m, "AbstractPhaseBasedCellCycleModel")
        .def(
            "ReadyToDivide",
            (bool(AbstractPhaseBasedCellCycleModel::*)()) &AbstractPhaseBasedCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "ResetForDivision",
            (void(AbstractPhaseBasedCellCycleModel::*)()) &AbstractPhaseBasedCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "UpdateCellCyclePhase",
            (void(AbstractPhaseBasedCellCycleModel::*)()) &AbstractPhaseBasedCellCycleModel::UpdateCellCyclePhase,
            " "  )
        .def(
            "GetCurrentCellCyclePhase",
            (::CellCyclePhase(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetCurrentCellCyclePhase,
            " "  )
        .def(
            "GetG1Duration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetG1Duration,
            " "  )
        .def(
            "GetStemCellG1Duration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetStemCellG1Duration,
            " "  )
        .def(
            "GetTransitCellG1Duration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetTransitCellG1Duration,
            " "  )
        .def(
            "GetSG2MDuration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetSG2MDuration,
            " "  )
        .def(
            "GetSDuration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetSDuration,
            " "  )
        .def(
            "GetG2Duration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetG2Duration,
            " "  )
        .def(
            "GetMDuration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetMDuration,
            " "  )
        .def(
            "SetStemCellG1Duration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetStemCellG1Duration,
            " " , py::arg("stemCellG1Duration") )
        .def(
            "SetTransitCellG1Duration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetTransitCellG1Duration,
            " " , py::arg("transitCellG1Duration") )
        .def(
            "SetSDuration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetSDuration,
            " " , py::arg("sDuration") )
        .def(
            "SetG2Duration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetG2Duration,
            " " , py::arg("g2Duration") )
        .def(
            "SetMDuration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetMDuration,
            " " , py::arg("mDuration") )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(AbstractPhaseBasedCellCycleModel::*)()) &AbstractPhaseBasedCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(AbstractPhaseBasedCellCycleModel::*)()) &AbstractPhaseBasedCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "GetMinimumGapDuration",
            (double(AbstractPhaseBasedCellCycleModel::*)() const ) &AbstractPhaseBasedCellCycleModel::GetMinimumGapDuration,
            " "  )
        .def(
            "SetMinimumGapDuration",
            (void(AbstractPhaseBasedCellCycleModel::*)(double)) &AbstractPhaseBasedCellCycleModel::SetMinimumGapDuration,
            " " , py::arg("minimumGapDuration") )
        .def(
            "OutputCellCycleModelParameters",
            (void(AbstractPhaseBasedCellCycleModel::*)(::out_stream &)) &AbstractPhaseBasedCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
