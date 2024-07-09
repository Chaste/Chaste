#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NullSrnModel.hpp"

#include "NullSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef NullSrnModel NullSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class NullSrnModel_Overrides : public NullSrnModel{
    public:
    using NullSrnModel::NullSrnModel;
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            NullSrnModel,
            SimulateToCurrentTime,
            );
    }
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            NullSrnModel,
            CreateSrnModel,
            );
    }

};
void register_NullSrnModel_class(py::module &m){
py::class_<NullSrnModel , NullSrnModel_Overrides , boost::shared_ptr<NullSrnModel >  , AbstractSrnModel  >(m, "NullSrnModel")
        .def(py::init< >())
        .def(
            "SimulateToCurrentTime",
            (void(NullSrnModel::*)()) &NullSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(NullSrnModel::*)()) &NullSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
    ;
}
