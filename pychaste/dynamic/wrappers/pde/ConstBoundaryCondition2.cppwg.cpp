#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ConstBoundaryCondition.hpp"

#include "ConstBoundaryCondition2.cppwg.hpp"

namespace py = pybind11;
typedef ConstBoundaryCondition<2 > ConstBoundaryCondition2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ConstBoundaryCondition2_Overrides : public ConstBoundaryCondition2{
    public:
    using ConstBoundaryCondition2::ConstBoundaryCondition;
    double GetValue(::ChastePoint<2> const & rX) const  override {
        PYBIND11_OVERRIDE(
            double,
            ConstBoundaryCondition2,
            GetValue,
                    rX);
    }

};
void register_ConstBoundaryCondition2_class(py::module &m){
py::class_<ConstBoundaryCondition2 , ConstBoundaryCondition2_Overrides , boost::shared_ptr<ConstBoundaryCondition2 >  , AbstractBoundaryCondition<2>  >(m, "ConstBoundaryCondition2")
        .def(py::init<double const >(), py::arg("value"))
        .def(
            "GetValue",
            (double(ConstBoundaryCondition2::*)(::ChastePoint<2> const &) const ) &ConstBoundaryCondition2::GetValue,
            " " , py::arg("rX") )
    ;
}
