#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"

#include "AbstractOdeBasedCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeBasedCellCycleModel AbstractOdeBasedCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeBasedCellCycleModel_Overrides : public AbstractOdeBasedCellCycleModel{
    public:
    using AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel;
    void SetBirthTime(double birthTime) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedCellCycleModel,
            SetBirthTime,
                    birthTime);
    }
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractOdeBasedCellCycleModel,
            ReadyToDivide,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedCellCycleModel,
            ResetForDivision,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeBasedCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AbstractOdeBasedCellCycleModel_class(py::module &m){
py::class_<AbstractOdeBasedCellCycleModel , AbstractOdeBasedCellCycleModel_Overrides , boost::shared_ptr<AbstractOdeBasedCellCycleModel >  , AbstractCellCycleModel  >(m, "AbstractOdeBasedCellCycleModel")
        .def(
            "GetOdeStopTime",
            (double(AbstractOdeBasedCellCycleModel::*)()) &AbstractOdeBasedCellCycleModel::GetOdeStopTime,
            " "  )
        .def(
            "SetBirthTime",
            (void(AbstractOdeBasedCellCycleModel::*)(double)) &AbstractOdeBasedCellCycleModel::SetBirthTime,
            " " , py::arg("birthTime") )
        .def(
            "ReadyToDivide",
            (bool(AbstractOdeBasedCellCycleModel::*)()) &AbstractOdeBasedCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "ResetForDivision",
            (void(AbstractOdeBasedCellCycleModel::*)()) &AbstractOdeBasedCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(AbstractOdeBasedCellCycleModel::*)(::out_stream &)) &AbstractOdeBasedCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
