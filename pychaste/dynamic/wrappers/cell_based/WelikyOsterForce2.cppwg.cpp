#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "WelikyOsterForce.hpp"

#include "WelikyOsterForce2.cppwg.hpp"

namespace py = pybind11;
typedef WelikyOsterForce<2 > WelikyOsterForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class WelikyOsterForce2_Overrides : public WelikyOsterForce2{
    public:
    using WelikyOsterForce2::WelikyOsterForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            WelikyOsterForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_WelikyOsterForce2_class(py::module &m){
py::class_<WelikyOsterForce2 , WelikyOsterForce2_Overrides , boost::shared_ptr<WelikyOsterForce2 >  , AbstractForce<2>  >(m, "WelikyOsterForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(WelikyOsterForce2::*)(::AbstractCellPopulation<2> &)) &WelikyOsterForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "GetWelikyOsterAreaParameter",
            (double(WelikyOsterForce2::*)()) &WelikyOsterForce2::GetWelikyOsterAreaParameter,
            " "  )
        .def(
            "GetWelikyOsterPerimeterParameter",
            (double(WelikyOsterForce2::*)()) &WelikyOsterForce2::GetWelikyOsterPerimeterParameter,
            " "  )
        .def(
            "SetWelikyOsterAreaParameter",
            (void(WelikyOsterForce2::*)(double)) &WelikyOsterForce2::SetWelikyOsterAreaParameter,
            " " , py::arg("welikyOsterAreaParameter") )
        .def(
            "SetWelikyOsterPerimeterParameter",
            (void(WelikyOsterForce2::*)(double)) &WelikyOsterForce2::SetWelikyOsterPerimeterParameter,
            " " , py::arg("welikyOsterPerimeterParameter") )
        .def(
            "OutputForceParameters",
            (void(WelikyOsterForce2::*)(::out_stream &)) &WelikyOsterForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
