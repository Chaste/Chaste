#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ContactInhibitionCellCycleModel.hpp"

#include "ContactInhibitionCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef ContactInhibitionCellCycleModel ContactInhibitionCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class ContactInhibitionCellCycleModel_Overrides : public ContactInhibitionCellCycleModel{
    public:
    using ContactInhibitionCellCycleModel::ContactInhibitionCellCycleModel;
    void UpdateCellCyclePhase() override {
        PYBIND11_OVERRIDE(
            void,
            ContactInhibitionCellCycleModel,
            UpdateCellCyclePhase,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            ContactInhibitionCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ContactInhibitionCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_ContactInhibitionCellCycleModel_class(py::module &m){
py::class_<ContactInhibitionCellCycleModel , ContactInhibitionCellCycleModel_Overrides , boost::shared_ptr<ContactInhibitionCellCycleModel >  , AbstractSimplePhaseBasedCellCycleModel  >(m, "ContactInhibitionCellCycleModel")
        .def(py::init< >())
        .def(
            "UpdateCellCyclePhase",
            (void(ContactInhibitionCellCycleModel::*)()) &ContactInhibitionCellCycleModel::UpdateCellCyclePhase,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(ContactInhibitionCellCycleModel::*)()) &ContactInhibitionCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "SetQuiescentVolumeFraction",
            (void(ContactInhibitionCellCycleModel::*)(double)) &ContactInhibitionCellCycleModel::SetQuiescentVolumeFraction,
            " " , py::arg("quiescentVolumeFraction") )
        .def(
            "GetQuiescentVolumeFraction",
            (double(ContactInhibitionCellCycleModel::*)() const ) &ContactInhibitionCellCycleModel::GetQuiescentVolumeFraction,
            " "  )
        .def(
            "SetEquilibriumVolume",
            (void(ContactInhibitionCellCycleModel::*)(double)) &ContactInhibitionCellCycleModel::SetEquilibriumVolume,
            " " , py::arg("equilibriumVolume") )
        .def(
            "GetEquilibriumVolume",
            (double(ContactInhibitionCellCycleModel::*)() const ) &ContactInhibitionCellCycleModel::GetEquilibriumVolume,
            " "  )
        .def(
            "GetCurrentQuiescentDuration",
            (double(ContactInhibitionCellCycleModel::*)() const ) &ContactInhibitionCellCycleModel::GetCurrentQuiescentDuration,
            " "  )
        .def(
            "GetCurrentQuiescentOnsetTime",
            (double(ContactInhibitionCellCycleModel::*)() const ) &ContactInhibitionCellCycleModel::GetCurrentQuiescentOnsetTime,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(ContactInhibitionCellCycleModel::*)(::out_stream &)) &ContactInhibitionCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
