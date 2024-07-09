#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AlwaysDivideCellCycleModel.hpp"

#include "AlwaysDivideCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AlwaysDivideCellCycleModel AlwaysDivideCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class AlwaysDivideCellCycleModel_Overrides : public AlwaysDivideCellCycleModel{
    public:
    using AlwaysDivideCellCycleModel::AlwaysDivideCellCycleModel;
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            AlwaysDivideCellCycleModel,
            ReadyToDivide,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            AlwaysDivideCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            AlwaysDivideCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            AlwaysDivideCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AlwaysDivideCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AlwaysDivideCellCycleModel_class(py::module &m){
py::class_<AlwaysDivideCellCycleModel , AlwaysDivideCellCycleModel_Overrides , boost::shared_ptr<AlwaysDivideCellCycleModel >  , AbstractCellCycleModel  >(m, "AlwaysDivideCellCycleModel")
        .def(py::init< >())
        .def(
            "ReadyToDivide",
            (bool(AlwaysDivideCellCycleModel::*)()) &AlwaysDivideCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(AlwaysDivideCellCycleModel::*)()) &AlwaysDivideCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "GetAverageTransitCellCycleTime",
            (double(AlwaysDivideCellCycleModel::*)()) &AlwaysDivideCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(AlwaysDivideCellCycleModel::*)()) &AlwaysDivideCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(AlwaysDivideCellCycleModel::*)(::out_stream &)) &AlwaysDivideCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
