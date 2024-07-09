#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BernoulliTrialCellCycleModel.hpp"

#include "BernoulliTrialCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef BernoulliTrialCellCycleModel BernoulliTrialCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class BernoulliTrialCellCycleModel_Overrides : public BernoulliTrialCellCycleModel{
    public:
    using BernoulliTrialCellCycleModel::BernoulliTrialCellCycleModel;
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            BernoulliTrialCellCycleModel,
            ReadyToDivide,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            BernoulliTrialCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            BernoulliTrialCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            BernoulliTrialCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BernoulliTrialCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_BernoulliTrialCellCycleModel_class(py::module &m){
py::class_<BernoulliTrialCellCycleModel , BernoulliTrialCellCycleModel_Overrides , boost::shared_ptr<BernoulliTrialCellCycleModel >  , AbstractCellCycleModel  >(m, "BernoulliTrialCellCycleModel")
        .def(py::init< >())
        .def(
            "ReadyToDivide",
            (bool(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "SetDivisionProbability",
            (void(BernoulliTrialCellCycleModel::*)(double)) &BernoulliTrialCellCycleModel::SetDivisionProbability,
            " " , py::arg("divisionProbability") )
        .def(
            "GetDivisionProbability",
            (double(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::GetDivisionProbability,
            " "  )
        .def(
            "SetMinimumDivisionAge",
            (void(BernoulliTrialCellCycleModel::*)(double)) &BernoulliTrialCellCycleModel::SetMinimumDivisionAge,
            " " , py::arg("minimumDivisionAge") )
        .def(
            "GetMinimumDivisionAge",
            (double(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::GetMinimumDivisionAge,
            " "  )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(BernoulliTrialCellCycleModel::*)()) &BernoulliTrialCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(BernoulliTrialCellCycleModel::*)(::out_stream &)) &BernoulliTrialCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
