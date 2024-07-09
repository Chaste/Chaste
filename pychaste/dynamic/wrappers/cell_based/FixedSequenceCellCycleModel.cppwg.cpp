#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedSequenceCellCycleModel.hpp"

#include "FixedSequenceCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef FixedSequenceCellCycleModel FixedSequenceCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class FixedSequenceCellCycleModel_Overrides : public FixedSequenceCellCycleModel{
    public:
    using FixedSequenceCellCycleModel::FixedSequenceCellCycleModel;
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            FixedSequenceCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void SetRate(double rate) override {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetRate,
                    rate);
    }
    double GetRate() override {
        PYBIND11_OVERRIDE(
            double,
            FixedSequenceCellCycleModel,
            GetRate,
            );
    }
    void SetStemCellG1Duration(double stemCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetStemCellG1Duration,
                    stemCellG1Duration);
    }
    void SetTransitCellG1Duration(double transitCellG1Duration) override {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetTransitCellG1Duration,
                    transitCellG1Duration);
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }
    void SetG1Duration() override {
        PYBIND11_OVERRIDE(
            void,
            FixedSequenceCellCycleModel,
            SetG1Duration,
            );
    }

};
void register_FixedSequenceCellCycleModel_class(py::module &m){
py::class_<FixedSequenceCellCycleModel , FixedSequenceCellCycleModel_Overrides , boost::shared_ptr<FixedSequenceCellCycleModel >  , ExponentialG1GenerationalCellCycleModel  >(m, "FixedSequenceCellCycleModel")
        .def(py::init< >())
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(FixedSequenceCellCycleModel::*)()) &FixedSequenceCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "SetRate",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetRate,
            " " , py::arg("rate") )
        .def(
            "GetRate",
            (double(FixedSequenceCellCycleModel::*)()) &FixedSequenceCellCycleModel::GetRate,
            " "  )
        .def(
            "SetStemCellG1Duration",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetStemCellG1Duration,
            " " , py::arg("stemCellG1Duration") )
        .def(
            "SetTransitCellG1Duration",
            (void(FixedSequenceCellCycleModel::*)(double)) &FixedSequenceCellCycleModel::SetTransitCellG1Duration,
            " " , py::arg("transitCellG1Duration") )
        .def(
            "OutputCellCycleModelParameters",
            (void(FixedSequenceCellCycleModel::*)(::out_stream &)) &FixedSequenceCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
