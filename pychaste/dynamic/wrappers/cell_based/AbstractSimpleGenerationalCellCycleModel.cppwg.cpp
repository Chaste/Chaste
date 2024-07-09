#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"

#include "AbstractSimpleGenerationalCellCycleModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractSimpleGenerationalCellCycleModel AbstractSimpleGenerationalCellCycleModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractSimpleGenerationalCellCycleModel_Overrides : public AbstractSimpleGenerationalCellCycleModel{
    public:
    using AbstractSimpleGenerationalCellCycleModel::AbstractSimpleGenerationalCellCycleModel;
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSimpleGenerationalCellCycleModel,
            ResetForDivision,
            );
    }
    void InitialiseDaughterCell() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSimpleGenerationalCellCycleModel,
            InitialiseDaughterCell,
            );
    }
    void OutputCellCycleModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractSimpleGenerationalCellCycleModel,
            OutputCellCycleModelParameters,
                    rParamsFile);
    }

};
void register_AbstractSimpleGenerationalCellCycleModel_class(py::module &m){
py::class_<AbstractSimpleGenerationalCellCycleModel , AbstractSimpleGenerationalCellCycleModel_Overrides , boost::shared_ptr<AbstractSimpleGenerationalCellCycleModel >  , AbstractSimplePhaseBasedCellCycleModel  >(m, "AbstractSimpleGenerationalCellCycleModel")
        .def(
            "ResetForDivision",
            (void(AbstractSimpleGenerationalCellCycleModel::*)()) &AbstractSimpleGenerationalCellCycleModel::ResetForDivision,
            " "  )
        .def(
            "InitialiseDaughterCell",
            (void(AbstractSimpleGenerationalCellCycleModel::*)()) &AbstractSimpleGenerationalCellCycleModel::InitialiseDaughterCell,
            " "  )
        .def(
            "SetGeneration",
            (void(AbstractSimpleGenerationalCellCycleModel::*)(unsigned int)) &AbstractSimpleGenerationalCellCycleModel::SetGeneration,
            " " , py::arg("generation") )
        .def(
            "GetGeneration",
            (unsigned int(AbstractSimpleGenerationalCellCycleModel::*)() const ) &AbstractSimpleGenerationalCellCycleModel::GetGeneration,
            " "  )
        .def(
            "SetMaxTransitGenerations",
            (void(AbstractSimpleGenerationalCellCycleModel::*)(unsigned int)) &AbstractSimpleGenerationalCellCycleModel::SetMaxTransitGenerations,
            " " , py::arg("maxTransitGenerations") )
        .def(
            "GetMaxTransitGenerations",
            (unsigned int(AbstractSimpleGenerationalCellCycleModel::*)() const ) &AbstractSimpleGenerationalCellCycleModel::GetMaxTransitGenerations,
            " "  )
        .def(
            "OutputCellCycleModelParameters",
            (void(AbstractSimpleGenerationalCellCycleModel::*)(::out_stream &)) &AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
