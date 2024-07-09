#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DiffusionForce.hpp"

#include "DiffusionForce2.cppwg.hpp"

namespace py = pybind11;
typedef DiffusionForce<2 > DiffusionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DiffusionForce2_Overrides : public DiffusionForce2{
    public:
    using DiffusionForce2::DiffusionForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            DiffusionForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            DiffusionForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_DiffusionForce2_class(py::module &m){
py::class_<DiffusionForce2 , DiffusionForce2_Overrides , boost::shared_ptr<DiffusionForce2 >  , AbstractForce<2>  >(m, "DiffusionForce2")
        .def(py::init< >())
        .def(
            "SetAbsoluteTemperature",
            (void(DiffusionForce2::*)(double)) &DiffusionForce2::SetAbsoluteTemperature,
            " " , py::arg("absoluteTemperature") )
        .def(
            "SetViscosity",
            (void(DiffusionForce2::*)(double)) &DiffusionForce2::SetViscosity,
            " " , py::arg("viscosity") )
        .def(
            "GetAbsoluteTemperature",
            (double(DiffusionForce2::*)()) &DiffusionForce2::GetAbsoluteTemperature,
            " "  )
        .def(
            "GetViscosity",
            (double(DiffusionForce2::*)()) &DiffusionForce2::GetViscosity,
            " "  )
        .def(
            "GetDiffusionScalingConstant",
            (double(DiffusionForce2::*)()) &DiffusionForce2::GetDiffusionScalingConstant,
            " "  )
        .def(
            "AddForceContribution",
            (void(DiffusionForce2::*)(::AbstractCellPopulation<2> &)) &DiffusionForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(DiffusionForce2::*)(::out_stream &)) &DiffusionForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
