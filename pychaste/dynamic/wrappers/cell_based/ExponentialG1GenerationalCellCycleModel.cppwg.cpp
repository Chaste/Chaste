#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"

#include "ExponentialG1GenerationalCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef ExponentialG1GenerationalCellCycleModel ExponentialG1GenerationalCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class ExponentialG1GenerationalCellCycleModel_Overrides : public ExponentialG1GenerationalCellCycleModel{
    public:
    using ExponentialG1GenerationalCellCycleModel::ExponentialG1GenerationalCellCycleModel;
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            ExponentialG1GenerationalCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetRate() override {
        PYBIND11_OVERRIDE(
            double,
            ExponentialG1GenerationalCellCycleModel,
            GetRate,
            );
    }
    void SetRate(double rate) override {
        PYBIND11_OVERRIDE(
            void,
            ExponentialG1GenerationalCellCycleModel,
            SetRate,
                    rate);
    }
    void SetStemCellG1Duration(double stemCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            ExponentialG1GenerationalCellCycleModel,
            SetStemCellG1Duration,
                    stemCellG1Duration);
    }
    void SetTransitCellG1Duration(double transitCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            ExponentialG1GenerationalCellCycleModel,
            SetTransitCellG1Duration,
                    transitCellG1Duration);
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ExponentialG1GenerationalCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }
    void SetG1Duration() override {
        PYBIND11_OVERRIDE(
            void,
            ExponentialG1GenerationalCellCycleModel,
            SetG1Duration,
            );
    }

};
void register_ExponentialG1GenerationalCellCycleModel_class(py::module &m){
py::class_<ExponentialG1GenerationalCellCycleModel , ExponentialG1GenerationalCellCycleModel_Overrides , boost::shared_ptr<ExponentialG1GenerationalCellCycleModel >  , AbstractSimpleGenerationalCellCycleModel  >(m, "ExponentialG1GenerationalCellCycleModel")
        .def(py::init< >())
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(ExponentialG1GenerationalCellCycleModel::*)()) &ExponentialG1GenerationalCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "GetRate",
            (double(ExponentialG1GenerationalCellCycleModel::*)()) &ExponentialG1GenerationalCellCycleModel::GetRate,
            " "  )
        .def(
            "SetRate",
            (void(ExponentialG1GenerationalCellCycleModel::*)(double)) &ExponentialG1GenerationalCellCycleModel::SetRate,
            " " , py::arg("rate") )
        .def(
            "SetStemCellG1Duration",
            (void(ExponentialG1GenerationalCellCycleModel::*)(double)) &ExponentialG1GenerationalCellCycleModel::SetStemCellG1Duration,
            " " , py::arg("stemCellG1Duration") )
        .def(
            "SetTransitCellG1Duration",
            (void(ExponentialG1GenerationalCellCycleModel::*)(double)) &ExponentialG1GenerationalCellCycleModel::SetTransitCellG1Duration,
            " " , py::arg("transitCellG1Duration") )
        .def(
            "OutputCellCycleModelParameters",
            (void(ExponentialG1GenerationalCellCycleModel::*)(::out_stream &)) &ExponentialG1GenerationalCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
