#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeCompressionForce.hpp"

#include "BuskeCompressionForce3.cppwg.hpp"

namespace py = pybind11;
typedef BuskeCompressionForce<3 > BuskeCompressionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class BuskeCompressionForce3_Overrides : public BuskeCompressionForce3{
    public:
    using BuskeCompressionForce3::BuskeCompressionForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeCompressionForce3_class(py::module &m){
py::class_<BuskeCompressionForce3 , BuskeCompressionForce3_Overrides , boost::shared_ptr<BuskeCompressionForce3 >  , AbstractForce<3>  >(m, "BuskeCompressionForce3")
        .def(py::init< >())
        .def(
            "GetCompressionEnergyParameter",
            (double(BuskeCompressionForce3::*)()) &BuskeCompressionForce3::GetCompressionEnergyParameter,
            " "  )
        .def(
            "SetCompressionEnergyParameter",
            (void(BuskeCompressionForce3::*)(double)) &BuskeCompressionForce3::SetCompressionEnergyParameter,
            " " , py::arg("compressionEnergyParameter") )
        .def(
            "AddForceContribution",
            (void(BuskeCompressionForce3::*)(::AbstractCellPopulation<3> &)) &BuskeCompressionForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(BuskeCompressionForce3::*)(::out_stream &)) &BuskeCompressionForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
