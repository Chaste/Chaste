#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellCycleModel.hpp"

#include "AbstractCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellCycleModel AbstractCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class AbstractCellCycleModel_Overrides : public AbstractCellCycleModel{
    public:
    using AbstractCellCycleModel::AbstractCellCycleModel;
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            Initialise,
            );
    }
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            InitialiseDaughterCell,
            );
    }
    void SetBirthTime(double birthTime) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            SetBirthTime,
                    birthTime);
    }
    bool ReadyToDivide() override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCellCycleModel,
            ReadyToDivide,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellCycleModel,
            ResetForDivision,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE_PURE(
            _AbstractCellCycleModelPtr,
            AbstractCellCycleModel,
            CreateCellCycleModel,
            );
    }
    bool CanCellTerminallyDifferentiate() override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCellCycleModel,
            CanCellTerminallyDifferentiate,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AbstractCellCycleModel_class(py::module &m){
py::class_<AbstractCellCycleModel , AbstractCellCycleModel_Overrides , boost::shared_ptr<AbstractCellCycleModel >   >(m, "AbstractCellCycleModel")
        .def(py::init< >())
        .def(
            "SetCell",
            (void(AbstractCellCycleModel::*)(::CellPtr)) &AbstractCellCycleModel::SetCell,
            " " , py::arg("pCell") )
        .def(
            "Initialise",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::Initialise,
            " "  )
        .def(
            "InitialiseDaughterCell",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::InitialiseDaughterCell,
            " "  )
        .def(
            "GetCell",
            (::CellPtr(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetCell,
            " "  )
        .def(
            "SetBirthTime",
            (void(AbstractCellCycleModel::*)(double)) &AbstractCellCycleModel::SetBirthTime,
            " " , py::arg("birthTime") )
        .def(
            "SetDimension",
            (void(AbstractCellCycleModel::*)(unsigned int)) &AbstractCellCycleModel::SetDimension,
            " " , py::arg("dimension") )
        .def(
            "GetDimension",
            (unsigned int(AbstractCellCycleModel::*)() const ) &AbstractCellCycleModel::GetDimension,
            " "  )
        .def(
            "GetBirthTime",
            (double(AbstractCellCycleModel::*)() const ) &AbstractCellCycleModel::GetBirthTime,
            " "  )
        .def(
            "GetAge",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAge,
            " "  )
        .def(
            "ReadyToDivide",
            (bool(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::ReadyToDivide,
            " "  )
        .def(
            "ResetForDivision",
            (void(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "CanCellTerminallyDifferentiate",
            (bool(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::CanCellTerminallyDifferentiate,
            " "  )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(AbstractCellCycleModel::*)()) &AbstractCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "OutputCellCycleModelInfo",
            (void(AbstractCellCycleModel::*)(::out_stream &)) &AbstractCellCycleModel::OutputCellCycleModelInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputCellCycleModelParameters",
            (void(AbstractCellCycleModel::*)(::out_stream &)) &AbstractCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
