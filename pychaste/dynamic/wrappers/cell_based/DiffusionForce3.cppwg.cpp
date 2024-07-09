#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DiffusionForce.hpp"

#include "DiffusionForce3.cppwg.hpp"

namespace py = pybind11;
typedef DiffusionForce<3 > DiffusionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DiffusionForce3_Overrides : public DiffusionForce3{
    public:
    using DiffusionForce3::DiffusionForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DiffusionForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DiffusionForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_DiffusionForce3_class(py::module &m){
py::class_<DiffusionForce3 , DiffusionForce3_Overrides , boost::shared_ptr<DiffusionForce3 >  , AbstractForce<3>  >(m, "DiffusionForce3")
        .def(py::init< >())
        .def(
            "SetAbsoluteTemperature",
            (void(DiffusionForce3::*)(double)) &DiffusionForce3::SetAbsoluteTemperature,
            " " , py::arg("absoluteTemperature") )
        .def(
            "SetViscosity",
            (void(DiffusionForce3::*)(double)) &DiffusionForce3::SetViscosity,
            " " , py::arg("viscosity") )
        .def(
            "GetAbsoluteTemperature",
            (double(DiffusionForce3::*)()) &DiffusionForce3::GetAbsoluteTemperature,
            " "  )
        .def(
            "GetViscosity",
            (double(DiffusionForce3::*)()) &DiffusionForce3::GetViscosity,
            " "  )
        .def(
            "GetDiffusionScalingConstant",
            (double(DiffusionForce3::*)()) &DiffusionForce3::GetDiffusionScalingConstant,
            " "  )
        .def(
            "AddForceContribution",
            (void(DiffusionForce3::*)(::AbstractCellPopulation<3> &)) &DiffusionForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(DiffusionForce3::*)(::out_stream &)) &DiffusionForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
