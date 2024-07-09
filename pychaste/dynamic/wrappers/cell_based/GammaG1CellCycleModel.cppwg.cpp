#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "GammaG1CellCycleModel.hpp"

#include "GammaG1CellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef GammaG1CellCycleModel GammaG1CellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class GammaG1CellCycleModel_Overrides : public GammaG1CellCycleModel{
    public:
    using GammaG1CellCycleModel::GammaG1CellCycleModel;
    void SetG1Duration() override {
        PYBIND11_OVERRIDE(
            void,
            GammaG1CellCycleModel,
            SetG1Duration,
            );
    }
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            GammaG1CellCycleModel,
            CreateCellCycleModel,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            GammaG1CellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_GammaG1CellCycleModel_class(py::module &m){
py::class_<GammaG1CellCycleModel , GammaG1CellCycleModel_Overrides , boost::shared_ptr<GammaG1CellCycleModel >  , AbstractSimplePhaseBasedCellCycleModel  >(m, "GammaG1CellCycleModel")
        .def(py::init< >())
        .def(
            "SetG1Duration",
            (void(GammaG1CellCycleModel::*)()) &GammaG1CellCycleModel::SetG1Duration,
            " "  )
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(GammaG1CellCycleModel::*)()) &GammaG1CellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "SetShape",
            (void(GammaG1CellCycleModel::*)(double)) &GammaG1CellCycleModel::SetShape,
            " " , py::arg("shape") )
        .def(
            "SetScale",
            (void(GammaG1CellCycleModel::*)(double)) &GammaG1CellCycleModel::SetScale,
            " " , py::arg("scale") )
        .def(
            "GetShape",
            (double(GammaG1CellCycleModel::*)() const ) &GammaG1CellCycleModel::GetShape,
            " "  )
        .def(
            "GetScale",
            (double(GammaG1CellCycleModel::*)() const ) &GammaG1CellCycleModel::GetScale,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(GammaG1CellCycleModel::*)(::out_stream &)) &GammaG1CellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
