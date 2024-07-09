#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "WelikyOsterForce.hpp"

#include "WelikyOsterForce3.cppwg.hpp"

namespace py = pybind11;
typedef WelikyOsterForce<3 > WelikyOsterForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class WelikyOsterForce3_Overrides : public WelikyOsterForce3{
    public:
    using WelikyOsterForce3::WelikyOsterForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_WelikyOsterForce3_class(py::module &m){
py::class_<WelikyOsterForce3 , WelikyOsterForce3_Overrides , boost::shared_ptr<WelikyOsterForce3 >  , AbstractForce<3>  >(m, "WelikyOsterForce3")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(WelikyOsterForce3::*)(::AbstractCellPopulation<3> &)) &WelikyOsterForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetWelikyOsterAreaParameter",
            (double(WelikyOsterForce3::*)()) &WelikyOsterForce3::GetWelikyOsterAreaParameter,
            " "  )
        .def(
            "GetWelikyOsterPerimeterParameter",
            (double(WelikyOsterForce3::*)()) &WelikyOsterForce3::GetWelikyOsterPerimeterParameter,
            " "  )
        .def(
            "SetWelikyOsterAreaParameter",
            (void(WelikyOsterForce3::*)(double)) &WelikyOsterForce3::SetWelikyOsterAreaParameter,
            " " , py::arg("welikyOsterAreaParameter") )
        .def(
            "SetWelikyOsterPerimeterParameter",
            (void(WelikyOsterForce3::*)(double)) &WelikyOsterForce3::SetWelikyOsterPerimeterParameter,
            " " , py::arg("welikyOsterPerimeterParameter") )
        .def(
            "OutputForceParameters",
            (void(WelikyOsterForce3::*)(::out_stream &)) &WelikyOsterForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
