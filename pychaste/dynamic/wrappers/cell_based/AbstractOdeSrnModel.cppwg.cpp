#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOdeSrnModel.hpp"

#include "AbstractOdeSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOdeSrnModel AbstractOdeSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractOdeSrnModel_Overrides : public AbstractOdeSrnModel{
    public:
    using AbstractOdeSrnModel::AbstractOdeSrnModel;
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeSrnModel,
            SimulateToCurrentTime,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeSrnModel,
            ResetForDivision,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeSrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }
    void ScaleSrnVariables(double const theta) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOdeSrnModel,
            ScaleSrnVariables,
                    theta);
    }

};
void register_AbstractOdeSrnModel_class(py::module &m){
py::class_<AbstractOdeSrnModel , AbstractOdeSrnModel_Overrides , boost::shared_ptr<AbstractOdeSrnModel >  , AbstractSrnModel  >(m, "AbstractOdeSrnModel")
        .def(
            "SimulateToCurrentTime",
            (void(AbstractOdeSrnModel::*)()) &AbstractOdeSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "ResetForDivision",
            (void(AbstractOdeSrnModel::*)()) &AbstractOdeSrnModel::ResetForDivision,
            " "  )
        .def(
            "SetInitialConditions",
            (void(AbstractOdeSrnModel::*)(::std::vector<double>)) &AbstractOdeSrnModel::SetInitialConditions,
            " " , py::arg("initialConditions") )
        .def(
            "OutputSrnModelParameters",
            (void(AbstractOdeSrnModel::*)(::out_stream &)) &AbstractOdeSrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "ScaleSrnVariables",
            (void(AbstractOdeSrnModel::*)(double const)) &AbstractOdeSrnModel::ScaleSrnVariables,
            " " , py::arg("theta") )
    ;
}
