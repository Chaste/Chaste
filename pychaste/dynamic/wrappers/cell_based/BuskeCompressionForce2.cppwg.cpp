#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeCompressionForce.hpp"

#include "BuskeCompressionForce2.cppwg.hpp"

namespace py = pybind11;
typedef BuskeCompressionForce<2 > BuskeCompressionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class BuskeCompressionForce2_Overrides : public BuskeCompressionForce2{
    public:
    using BuskeCompressionForce2::BuskeCompressionForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeCompressionForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeCompressionForce2_class(py::module &m){
py::class_<BuskeCompressionForce2 , BuskeCompressionForce2_Overrides , boost::shared_ptr<BuskeCompressionForce2 >  , AbstractForce<2>  >(m, "BuskeCompressionForce2")
        .def(py::init< >())
        .def(
            "GetCompressionEnergyParameter",
            (double(BuskeCompressionForce2::*)()) &BuskeCompressionForce2::GetCompressionEnergyParameter,
            " "  )
        .def(
            "SetCompressionEnergyParameter",
            (void(BuskeCompressionForce2::*)(double)) &BuskeCompressionForce2::SetCompressionEnergyParameter,
            " " , py::arg("compressionEnergyParameter") )
        .def(
            "AddForceContribution",
            (void(BuskeCompressionForce2::*)(::AbstractCellPopulation<2> &)) &BuskeCompressionForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(BuskeCompressionForce2::*)(::out_stream &)) &BuskeCompressionForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
