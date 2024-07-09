#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchSrnModel.hpp"

#include "DeltaNotchSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchSrnModel DeltaNotchSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class DeltaNotchSrnModel_Overrides : public DeltaNotchSrnModel{
    public:
    using DeltaNotchSrnModel::DeltaNotchSrnModel;
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            DeltaNotchSrnModel,
            CreateSrnModel,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchSrnModel,
            Initialise,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchSrnModel,
            SimulateToCurrentTime,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchSrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }

};
void register_DeltaNotchSrnModel_class(py::module &m){
py::class_<DeltaNotchSrnModel , DeltaNotchSrnModel_Overrides , boost::shared_ptr<DeltaNotchSrnModel >  , AbstractOdeSrnModel  >(m, "DeltaNotchSrnModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "Initialise",
            (void(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::Initialise,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "UpdateDeltaNotch",
            (void(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::UpdateDeltaNotch,
            " "  )
        .def(
            "GetNotch",
            (double(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::GetNotch,
            " "  )
        .def(
            "GetDelta",
            (double(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::GetDelta,
            " "  )
        .def(
            "GetMeanNeighbouringDelta",
            (double(DeltaNotchSrnModel::*)()) &DeltaNotchSrnModel::GetMeanNeighbouringDelta,
            " "  )
        .def(
            "OutputSrnModelParameters",
            (void(DeltaNotchSrnModel::*)(::out_stream &)) &DeltaNotchSrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
    ;
}
