#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "TysonNovakCellCycleModel.hpp"

#include "TysonNovakCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef TysonNovakCellCycleModel TysonNovakCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class TysonNovakCellCycleModel_Overrides : public TysonNovakCellCycleModel{
    public:
    using TysonNovakCellCycleModel::TysonNovakCellCycleModel;
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            TysonNovakCellCycleModel,
            Initialise,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            TysonNovakCellCycleModel,
            ResetForDivision,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            TysonNovakCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            TysonNovakCellCycleModel,
            InitialiseDaughterCell,
            );
    }
    double GetAverageTransitCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            TysonNovakCellCycleModel,
            GetAverageTransitCellCycleTime,
            );
    }
    double GetAverageStemCellCycleTime() override {
        PYBIND11_OVERRIDE(
            double,
            TysonNovakCellCycleModel,
            GetAverageStemCellCycleTime,
            );
    }
    bool CanCellTerminallyDifferentiate() override {
        PYBIND11_OVERRIDE(
            bool,
            TysonNovakCellCycleModel,
            CanCellTerminallyDifferentiate,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            TysonNovakCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_TysonNovakCellCycleModel_class(py::module &m){
py::class_<TysonNovakCellCycleModel , TysonNovakCellCycleModel_Overrides , boost::shared_ptr<TysonNovakCellCycleModel >  , AbstractOdeBasedCellCycleModel  >(m, "TysonNovakCellCycleModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "Initialise",
            (void(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::Initialise,
            " "  )
        .def(
            "ResetForDivision",
            (void(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "InitialiseDaughterCell",
            (void(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::InitialiseDaughterCell,
            " "  )
        .def(
            "GetAverageTransitCellCycleTime",
            (double(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::GetAverageTransitCellCycleTime,
            " "  )
        .def(
            "GetAverageStemCellCycleTime",
            (double(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::GetAverageStemCellCycleTime,
            " "  )
        .def(
            "CanCellTerminallyDifferentiate",
            (bool(TysonNovakCellCycleModel::*)()) &TysonNovakCellCycleModel::CanCellTerminallyDifferentiate,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(TysonNovakCellCycleModel::*)(::out_stream &)) &TysonNovakCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
