#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NoCellCycleModel.hpp"

#include "NoCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef NoCellCycleModel NoCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class NoCellCycleModel_Overrides : public NoCellCycleModel{
    public:
    using NoCellCycleModel::NoCellCycleModel;
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            NoCellCycleModel,
            ReadyToDivide,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            NoCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            NoCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            NoCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            NoCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_NoCellCycleModel_class(py::module &m){
py::class_<NoCellCycleModel , NoCellCycleModel_Overrides , boost::shared_ptr<NoCellCycleModel >  , AbstractCellCycleModel  >(m, "NoCellCycleModel")
        .def(py::init< >())
        .def(
            "ReadyToDivide",
            (bool(NoCellCycleModel::*)()) &NoCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(NoCellCycleModel::*)()) &NoCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "GetAverageTransitCellCycleTime",
            (double(NoCellCycleModel::*)()) &NoCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(NoCellCycleModel::*)()) &NoCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(NoCellCycleModel::*)(::out_stream &)) &NoCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
