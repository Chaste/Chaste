#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"

#include "Alarcon2004OxygenBasedCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef Alarcon2004OxygenBasedCellCycleModel Alarcon2004OxygenBasedCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class Alarcon2004OxygenBasedCellCycleModel_Overrides : public Alarcon2004OxygenBasedCellCycleModel{
    public:
    using Alarcon2004OxygenBasedCellCycleModel::Alarcon2004OxygenBasedCellCycleModel;
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            Alarcon2004OxygenBasedCellCycleModel,
            ResetForDivision,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            Alarcon2004OxygenBasedCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            Alarcon2004OxygenBasedCellCycleModel,
            Initialise,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            Alarcon2004OxygenBasedCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_Alarcon2004OxygenBasedCellCycleModel_class(py::module &m){
py::class_<Alarcon2004OxygenBasedCellCycleModel , Alarcon2004OxygenBasedCellCycleModel_Overrides , boost::shared_ptr<Alarcon2004OxygenBasedCellCycleModel >  , AbstractOdeBasedPhaseBasedCellCycleModel  >(m, "Alarcon2004OxygenBasedCellCycleModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "ResetForDivision",
            (void(Alarcon2004OxygenBasedCellCycleModel::*)()) &Alarcon2004OxygenBasedCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(Alarcon2004OxygenBasedCellCycleModel::*)()) &Alarcon2004OxygenBasedCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "Initialise",
            (void(Alarcon2004OxygenBasedCellCycleModel::*)()) &Alarcon2004OxygenBasedCellCycleModel::Initialise,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(Alarcon2004OxygenBasedCellCycleModel::*)(::out_stream &)) &Alarcon2004OxygenBasedCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
