#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "UniformCellCycleModel.hpp"

#include "UniformCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef UniformCellCycleModel UniformCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class UniformCellCycleModel_Overrides : public UniformCellCycleModel{
    public:
    using UniformCellCycleModel::UniformCellCycleModel;
    void SetCellCycleDuration() override {
        PYBIND11_OVERRIDE(
            void,
            UniformCellCycleModel,
            SetCellCycleDuration,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            UniformCellCycleModel,
            CreateCellCycleModel,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            UniformCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            UniformCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            UniformCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_UniformCellCycleModel_class(py::module &m){
py::class_<UniformCellCycleModel , UniformCellCycleModel_Overrides , boost::shared_ptr<UniformCellCycleModel >  , AbstractSimpleCellCycleModel  >(m, "UniformCellCycleModel")
        .def(py::init< >())
        .def(
            "SetCellCycleDuration",
            (void(UniformCellCycleModel::*)()) &UniformCellCycleModel::SetCellCycleDuration,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(UniformCellCycleModel::*)()) &UniformCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "GetMinCellCycleDuration",
            (double(UniformCellCycleModel::*)()) &UniformCellCycleModel::GetMinCellCycleDuration,
            " "  )
        .def(
            "SetMinCellCycleDuration",
            (void(UniformCellCycleModel::*)(double)) &UniformCellCycleModel::SetMinCellCycleDuration,
            " " , py::arg("minCellCycleDuration") )
        .def(
            "GetMaxCellCycleDuration",
            (double(UniformCellCycleModel::*)()) &UniformCellCycleModel::GetMaxCellCycleDuration,
            " "  )
        .def(
            "SetMaxCellCycleDuration",
            (void(UniformCellCycleModel::*)(double)) &UniformCellCycleModel::SetMaxCellCycleDuration,
            " " , py::arg("maxCellCycleDuration") )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(UniformCellCycleModel::*)()) &UniformCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(UniformCellCycleModel::*)()) &UniformCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(UniformCellCycleModel::*)(::out_stream &)) &UniformCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
