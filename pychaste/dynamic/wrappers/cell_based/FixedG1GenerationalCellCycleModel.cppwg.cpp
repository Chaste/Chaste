#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "FixedG1GenerationalCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef FixedG1GenerationalCellCycleModel FixedG1GenerationalCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractCellCycleModel * _AbstractCellCycleModelPtr;

class FixedG1GenerationalCellCycleModel_Overrides : public FixedG1GenerationalCellCycleModel{
    public:
    using FixedG1GenerationalCellCycleModel::FixedG1GenerationalCellCycleModel;
    ::AbstractCellCycleModel * CreateCellCycleModel() override {
        PYBIND11_OVERRIDE(
            _AbstractCellCycleModelPtr,
            FixedG1GenerationalCellCycleModel,
            CreateCellCycleModel,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            FixedG1GenerationalCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_FixedG1GenerationalCellCycleModel_class(py::module &m){
py::class_<FixedG1GenerationalCellCycleModel , FixedG1GenerationalCellCycleModel_Overrides , boost::shared_ptr<FixedG1GenerationalCellCycleModel >  , AbstractSimpleGenerationalCellCycleModel  >(m, "FixedG1GenerationalCellCycleModel")
        .def(py::init< >())
        .def(
            "CreateCellCycleModel",
            (::AbstractCellCycleModel *(FixedG1GenerationalCellCycleModel::*)()) &FixedG1GenerationalCellCycleModel::CreateCellCycleModel,
            " "  , py::return_value_policy::reference)
        .def(
            "OutputCellCycleModelParameters",
            (void(FixedG1GenerationalCellCycleModel::*)(::out_stream &)) &FixedG1GenerationalCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
