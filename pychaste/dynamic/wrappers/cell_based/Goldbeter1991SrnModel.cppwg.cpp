#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Goldbeter1991SrnModel.hpp"

#include "Goldbeter1991SrnModel.cppwg.hpp"

namespace py = pybind11;
typedef Goldbeter1991SrnModel Goldbeter1991SrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class Goldbeter1991SrnModel_Overrides : public Goldbeter1991SrnModel{
    public:
    using Goldbeter1991SrnModel::Goldbeter1991SrnModel;
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            Goldbeter1991SrnModel,
            CreateSrnModel,
            );
    }
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            Goldbeter1991SrnModel,
            Initialise,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            Goldbeter1991SrnModel,
            SimulateToCurrentTime,
            );
    }
    void OutputSrnModelParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            Goldbeter1991SrnModel,
            OutputSrnModelParameters,
                    rParamsFile);
    }

};
void register_Goldbeter1991SrnModel_class(py::module &m){
py::class_<Goldbeter1991SrnModel , Goldbeter1991SrnModel_Overrides , boost::shared_ptr<Goldbeter1991SrnModel >  , AbstractOdeSrnModel  >(m, "Goldbeter1991SrnModel")
        .def(py::init<::boost::shared_ptr<AbstractCellCycleModelOdeSolver> >(), py::arg("pOdeSolver") = boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "Initialise",
            (void(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::Initialise,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "OutputSrnModelParameters",
            (void(Goldbeter1991SrnModel::*)(::out_stream &)) &Goldbeter1991SrnModel::OutputSrnModelParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetC",
            (double(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::GetC,
            " "  )
        .def(
            "GetM",
            (double(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::GetM,
            " "  )
        .def(
            "GetX",
            (double(Goldbeter1991SrnModel::*)()) &Goldbeter1991SrnModel::GetX,
            " "  )
    ;
}
