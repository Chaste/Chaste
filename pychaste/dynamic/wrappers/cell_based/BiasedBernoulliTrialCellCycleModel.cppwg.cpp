#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BiasedBernoulliTrialCellCycleModel.hpp"

#include "BiasedBernoulliTrialCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef BiasedBernoulliTrialCellCycleModel BiasedBernoulliTrialCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class BiasedBernoulliTrialCellCycleModel_Overrides : public BiasedBernoulliTrialCellCycleModel{
    public:
    using BiasedBernoulliTrialCellCycleModel::BiasedBernoulliTrialCellCycleModel;
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            BiasedBernoulliTrialCellCycleModel,
            ReadyToDivide,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            BiasedBernoulliTrialCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            BiasedBernoulliTrialCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            BiasedBernoulliTrialCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BiasedBernoulliTrialCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_BiasedBernoulliTrialCellCycleModel_class(py::module &m){
py::class_<BiasedBernoulliTrialCellCycleModel , BiasedBernoulliTrialCellCycleModel_Overrides , boost::shared_ptr<BiasedBernoulliTrialCellCycleModel >  , AbstractCellCycleModel  >(m, "BiasedBernoulliTrialCellCycleModel")
        .def(py::init< >())
        .def(
            "ReadyToDivide",
            (bool(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "SetMaxDivisionProbability",
            (void(BiasedBernoulliTrialCellCycleModel::*)(double)) &BiasedBernoulliTrialCellCycleModel::SetMaxDivisionProbability,
            " " , py::arg("maxDivisionProbability") )
        .def(
            "GetMaxDivisionProbability",
            (double(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::GetMaxDivisionProbability,
            " "  )
        .def(
            "SetMinimumDivisionAge",
            (void(BiasedBernoulliTrialCellCycleModel::*)(double)) &BiasedBernoulliTrialCellCycleModel::SetMinimumDivisionAge,
            " " , py::arg("minimumDivisionAge") )
        .def(
            "GetMinimumDivisionAge",
            (double(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::GetMinimumDivisionAge,
            " "  )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(BiasedBernoulliTrialCellCycleModel::*)()) &BiasedBernoulliTrialCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(BiasedBernoulliTrialCellCycleModel::*)(::out_stream &)) &BiasedBernoulliTrialCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
