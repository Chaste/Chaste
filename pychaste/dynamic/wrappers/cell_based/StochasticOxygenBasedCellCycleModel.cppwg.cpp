#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"

#include "StochasticOxygenBasedCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef StochasticOxygenBasedCellCycleModel StochasticOxygenBasedCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class StochasticOxygenBasedCellCycleModel_Overrides : public StochasticOxygenBasedCellCycleModel{
    public:
    using StochasticOxygenBasedCellCycleModel::StochasticOxygenBasedCellCycleModel;
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            StochasticOxygenBasedCellCycleModel,
            InitialiseDaughterCell,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            StochasticOxygenBasedCellCycleModel,
            Initialise,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            StochasticOxygenBasedCellCycleModel,
            ResetForDivision,
            );
    }
    double GetG2Duration() const  override {
        PYBIND11_OVERRIDE(
            double,
            StochasticOxygenBasedCellCycleModel,
            GetG2Duration,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            StochasticOxygenBasedCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            StochasticOxygenBasedCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_StochasticOxygenBasedCellCycleModel_class(py::module &m){
py::class_<StochasticOxygenBasedCellCycleModel , StochasticOxygenBasedCellCycleModel_Overrides , boost::shared_ptr<StochasticOxygenBasedCellCycleModel >  , SimpleOxygenBasedCellCycleModel  >(m, "StochasticOxygenBasedCellCycleModel")
        .def(py::init< >())
        .def(
            "InitialiseDaughterCell",
            (void(StochasticOxygenBasedCellCycleModel::*)()) &StochasticOxygenBasedCellCycleModel::InitialiseDaughterCell,
            " "  )
        .def(
            "Initialise",
            (void(StochasticOxygenBasedCellCycleModel::*)()) &StochasticOxygenBasedCellCycleModel::Initialise,
            " "  )
        .def(
            "ResetForDivision",
            (void(StochasticOxygenBasedCellCycleModel::*)()) &StochasticOxygenBasedCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "GetG2Duration",
            (double(StochasticOxygenBasedCellCycleModel::*)() const ) &StochasticOxygenBasedCellCycleModel::GetG2Duration,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(StochasticOxygenBasedCellCycleModel::*)()) &StochasticOxygenBasedCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "OutputCellCycleModelParameters",
            (void(StochasticOxygenBasedCellCycleModel::*)(::out_stream &)) &StochasticOxygenBasedCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
