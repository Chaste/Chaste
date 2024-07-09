#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ConstBoundaryCondition.hpp"

#include "ConstBoundaryCondition3.cppwg.hpp"

namespace py = pybind11;
typedef ConstBoundaryCondition<3 > ConstBoundaryCondition3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ConstBoundaryCondition3_Overrides : public ConstBoundaryCondition3{
    public:
    using ConstBoundaryCondition3::ConstBoundaryCondition;
    double GetValue(::ChastePoint<3> const & rX) const  override {
        PYBIND11_OVERRIDE(
            double,
            ConstBoundaryCondition3,
            GetValue,
                    rX);
    }

};
void register_ConstBoundaryCondition3_class(py::module &m){
py::class_<ConstBoundaryCondition3 , ConstBoundaryCondition3_Overrides , boost::shared_ptr<ConstBoundaryCondition3 >  , AbstractBoundaryCondition<3>  >(m, "ConstBoundaryCondition3")
        .def(py::init<double const >(), py::arg("value"))
        .def(
            "GetValue",
            (double(ConstBoundaryCondition3::*)(::ChastePoint<3> const &) const ) &ConstBoundaryCondition3::GetValue,
            " " , py::arg("rX") )
    ;
}
